from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import sep
from scipy import ndimage

from .astrometry import estimate_pixel_scale_arcsec
from .calibration import calibrate_science_frame, load_master_dark, load_master_frame
from .detection import build_static_sky_model, create_residual, extract_residual_detections, robust_location_scale
from .fits_io import read_fits_image
from .metadata import inspect_sequence
from .registration import register_frame


@dataclass(frozen=True)
class DetectionOptimizationResult:
    settings: dict[str, float | int]
    metrics: dict[str, float | int | str]
    warnings: tuple[str, ...] = ()

    def summary(self) -> str:
        scale = self.metrics.get("pixel_scale_arcsec", "WCS yok")
        scale_text = f"{float(scale):.2f} arcsec/px" if isinstance(scale, (int, float)) else str(scale)
        return (
            f"{self.settings['sigma']:.1f} sigma · FWHM {self.metrics['fwhm_px']:.2f} px · "
            f"{scale_text} · {self.settings['min_frames']} kare"
        )


def optimize_detection_settings(
    fits_files: list[str],
    *,
    master_paths: dict[str, str] | None = None,
    progress_callback: Callable[[int, str], None] | None = None,
    stop_callback: Callable[[], bool] | None = None,
) -> DetectionOptimizationResult:
    if len(fits_files) < 3:
        raise ValueError("Otomatik optimizasyon için en az üç FITS karesi gerekir.")

    def emit(value: int, message: str):
        if progress_callback:
            progress_callback(int(value), message)
        if stop_callback and stop_callback():
            raise RuntimeError("Optimizasyon durduruldu.")

    emit(2, "Sekans metadatası ölçülüyor")
    frames = inspect_sequence(fits_files)
    if len(frames) < 3:
        raise ValueError("Analize uygun zaman sıralı en az üç kare bulunamadı.")
    if any("shape_mismatch" in frame.quality_flags for frame in frames):
        raise ValueError("Kare boyutları eşit olmadığı için otomatik profil üretilemedi.")

    masters, dark_exposure = _load_masters(master_paths or {})
    sample_count = min(7, len(frames))
    sample_indices = sorted(set(np.linspace(0, len(frames) - 1, sample_count, dtype=int).tolist()))
    sample_images: list[np.ndarray] = []
    headers = []
    noise_values = []
    fwhm_values = []
    density_values = []
    correlation_values = []

    for order, frame_index in enumerate(sample_indices, start=1):
        emit(5 + int(25 * order / len(sample_indices)), f"Görüntü kalitesi ölçülüyor {order}/{len(sample_indices)}")
        frame = frames[frame_index]
        raw, header = read_fits_image(frame.file_path, header=True)
        calibrated = calibrate_science_frame(
            raw,
            master_bias=masters["bias"],
            master_dark=masters["dark"],
            master_flat=masters["flat"],
            science_exposure=frame.exposure_seconds,
            dark_exposure=dark_exposure,
        )
        crop = _center_crop(calibrated, max_dim=1100)
        _, noise = robust_location_scale(crop)
        fwhm, density, correlation = _measure_frame_quality(crop, noise)
        sample_images.append(crop)
        headers.append(header)
        noise_values.append(noise)
        if fwhm is not None:
            fwhm_values.append(fwhm)
        density_values.append(density)
        correlation_values.append(correlation)

    fwhm_px = float(np.median(fwhm_values)) if fwhm_values else 3.0
    background_rms = float(np.median(noise_values))
    source_density = float(np.median(density_values))
    noise_correlation = float(np.median(correlation_values))
    pixel_scales = [estimate_pixel_scale_arcsec(header) for header in headers]
    pixel_scales = [value for value in pixel_scales if value is not None and np.isfinite(value) and value > 0]
    camera = frames[0].camera
    pixel_scale = float(np.median(pixel_scales)) if pixel_scales else camera.image_scale_arcsec_px

    emit(36, "Örnek kareler alt piksel hizalanıyor")
    reference_position = len(sample_images) // 2
    reference = sample_images[reference_position]
    aligned = []
    registration_rms = []
    for order, image in enumerate(sample_images, start=1):
        if order - 1 == reference_position:
            aligned.append(reference)
            registration_rms.append(0.0)
        else:
            registered, solution = register_frame(reference, image)
            aligned.append(registered)
            registration_rms.append(float(solution.rms_px))
        emit(36 + int(22 * order / len(sample_images)), f"Hizalama ölçülüyor {order}/{len(sample_images)}")

    median_registration_rms = float(np.median(registration_rms))
    static_model = build_static_sky_model(aligned)
    residuals = [create_residual(image, static_model) for image in aligned]
    min_pixels = int(np.clip(round(0.75 * np.pi * (fwhm_px / 2.355) ** 2), 5, 18))
    edge_margin = int(np.clip(round(max(8.0, 2.5 * fwhm_px)), 8, 32))
    sigma_candidates = [4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    candidate_limit = 8 if source_density < 700 else 12 if source_density < 2200 else 16
    selected_sigma = sigma_candidates[-1]
    selected_counts: list[int] = []
    scan_results: dict[float, tuple[list[int], float]] = {}

    for sigma_index, sigma in enumerate(sigma_candidates, start=1):
        counts = []
        artifact_count = 0
        detection_count = 0
        for residual in residuals:
            detections, _rms = extract_residual_detections(
                residual,
                0,
                threshold_sigma=sigma,
                min_pixels=min_pixels,
                max_sources=96,
                edge_margin=edge_margin,
                expected_fwhm_px=fwhm_px,
            )
            counts.append(len(detections))
            detection_count += len(detections)
            artifact_count += sum(bool(set(item.flags) & {"star_subtraction_dipole", "near_saturation", "hot_pixel_or_cosmic_ray"}) for item in detections)
        artifact_fraction = artifact_count / max(1, detection_count)
        scan_results[sigma] = (counts, artifact_fraction)
        p90 = float(np.percentile(counts, 90)) if counts else 0.0
        if p90 <= candidate_limit and artifact_fraction <= 0.45:
            selected_sigma = sigma
            selected_counts = counts
            break
        emit(60 + int(20 * sigma_index / len(sigma_candidates)), f"Residual eşiği deneniyor: {sigma:.1f} sigma")
    if not selected_counts:
        selected_counts = scan_results[selected_sigma][0]

    cadence_minutes = _median_cadence_minutes(frames)
    exposure_seconds = _median_positive([frame.exposure_seconds for frame in frames]) or 60.0
    min_motion_px, max_step_px = _motion_limits(
        pixel_scale_arcsec=pixel_scale,
        cadence_minutes=cadence_minutes,
        exposure_seconds=exposure_seconds,
        fwhm_px=fwhm_px,
        registration_rms_px=median_registration_rms,
    )
    frame_count = len(frames)
    min_frames = 3 if frame_count <= 5 else 4 if frame_count <= 15 else 5 if frame_count <= 40 else 6
    max_fit_rms = float(np.clip(0.25 * fwhm_px + 0.35 * median_registration_rms, 0.70, 1.60))
    match_tolerance = float(np.clip(0.65 * fwhm_px + 0.35 * median_registration_rms, 1.6, 3.5))
    p90_count = float(np.percentile(selected_counts, 90)) if selected_counts else 0.0
    max_sources = int(np.clip(np.ceil(max(candidate_limit, p90_count) * 1.35), 12, 32))
    min_occupancy = 0.50 if frame_count <= 7 else 0.58 if frame_count <= 20 else 0.65
    max_missing_gap = 2 if frame_count <= 20 else int(np.clip(round(frame_count * 0.07), 2, 4))

    warnings = []
    if pixel_scale is None:
        warnings.append("WCS piksel ölçeği bulunamadı; hareket sınırları PSF ve kadanstan türetildi.")
    if median_registration_rms > 1.5:
        warnings.append("Örnek hizalama RMS yüksek; plate solve ve takip kalitesi kontrol edilmeli.")
    if source_density > 2500:
        warnings.append("Alan kalabalık; residual eşik ve minimum kare sayısı yükseltildi.")
    if noise_correlation > 1.35:
        warnings.append("Korelasyonlu arka plan yapısı ölçüldü; agresif düşük sigma kullanılmadı.")
    if fwhm_px < 1.8:
        warnings.append("Görüntü örneklemesi düşük; tek piksel ve sıcak piksel filtresi sıkılaştırıldı.")

    settings: dict[str, float | int] = {
        "sigma": float(selected_sigma),
        "min_pixels": min_pixels,
        "min_frames": min_frames,
        "max_sources": max_sources,
        "edge_margin": edge_margin,
        "expected_fwhm_px": round(fwhm_px, 3),
        "min_motion_px_per_frame": round(min_motion_px, 4),
        "max_step_px": round(max_step_px, 3),
        "min_median_snr": round(max(5.0, selected_sigma * 1.02), 3),
        "max_fit_rms_px": round(max_fit_rms, 3),
        "strong_fit_rms_px": round(max(0.55, max_fit_rms * 0.62), 3),
        "match_tolerance_px": round(match_tolerance, 3),
        "min_track_occupancy": round(min_occupancy, 3),
        "max_missing_gap_frames": max_missing_gap,
        "max_artifact_fraction": 0.34,
        "persistence_fraction": round(float(np.clip(min_frames / max(frame_count, 1) * 1.4, 0.15, 0.35)), 3),
    }
    metrics: dict[str, float | int | str] = {
        "sample_frames": len(sample_images),
        "frame_count": frame_count,
        "fwhm_px": round(fwhm_px, 3),
        "background_rms": round(background_rms, 5),
        "source_density_mpix": round(source_density, 1),
        "noise_correlation": round(noise_correlation, 3),
        "registration_rms_px": round(median_registration_rms, 3),
        "pixel_scale_arcsec": round(pixel_scale, 4) if pixel_scale is not None else "WCS yok",
        "cadence_minutes": round(cadence_minutes, 4),
        "exposure_seconds": round(exposure_seconds, 3),
        "residuals_p90": round(p90_count, 1),
        "camera": camera.instrument or "header missing",
        "detector": camera.detector or "header missing",
        "pixel_size_um": round(camera.pixel_size_x_um, 3) if camera.pixel_size_x_um is not None else "header missing",
        "binning": f"{camera.binning_x}x{camera.binning_y}" if camera.binning_x and camera.binning_y else "header missing",
        "gain": round(camera.gain, 3) if camera.gain is not None else "header missing",
        "sensor_temperature_c": round(camera.sensor_temperature_c, 2) if camera.sensor_temperature_c is not None else "header missing",
        "focal_length_mm": round(camera.focal_length_mm, 2) if camera.focal_length_mm is not None else "header missing",
    }
    emit(100, "Otomatik tespit profili hazır")
    return DetectionOptimizationResult(settings=settings, metrics=metrics, warnings=tuple(warnings))


def _load_masters(paths: dict[str, str]):
    masters = {"bias": None, "dark": None, "flat": None}
    dark_exposure = None
    if paths.get("bias"):
        masters["bias"] = load_master_frame(paths["bias"])
    if paths.get("dark"):
        masters["dark"], dark_exposure = load_master_dark(paths["dark"])
    if paths.get("flat"):
        masters["flat"] = load_master_frame(paths["flat"])
    return masters, dark_exposure


def _center_crop(image: np.ndarray, *, max_dim: int) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    height, width = data.shape
    crop_height = min(height, max_dim)
    crop_width = min(width, max_dim)
    y0 = max(0, (height - crop_height) // 2)
    x0 = max(0, (width - crop_width) // 2)
    return np.ascontiguousarray(data[y0 : y0 + crop_height, x0 : x0 + crop_width])


def _measure_frame_quality(image: np.ndarray, noise: float):
    data = np.ascontiguousarray(image, dtype=np.float32)
    background = sep.Background(data)
    centered = data - background
    objects = sep.extract(centered, 6.0, err=max(float(background.globalrms), noise, 1e-6), minarea=3)
    height, width = data.shape
    fwhm = None
    valid_count = 0
    if len(objects):
        major = np.asarray(objects["a"], dtype=float)
        minor = np.asarray(objects["b"], dtype=float)
        valid = (
            (objects["x"] > 12)
            & (objects["x"] < width - 12)
            & (objects["y"] > 12)
            & (objects["y"] < height - 12)
            & (major > 0.55)
            & (minor > 0.55)
            & (major / np.maximum(minor, 1e-6) < 2.2)
        )
        widths = 2.355 * np.sqrt(major[valid] * minor[valid])
        widths = widths[np.isfinite(widths) & (widths >= 1.2) & (widths <= 15.0)]
        valid_count = int(np.count_nonzero(valid))
        if widths.size:
            fwhm = float(np.median(widths))
    density = valid_count / max(1e-6, data.size / 1_000_000.0)
    high_pass = centered - ndimage.gaussian_filter(centered, sigma=1.0)
    _hp_median, hp_noise = robust_location_scale(high_pass)
    correlation = float(np.clip(noise / max(hp_noise, 1e-6), 0.5, 3.0))
    return fwhm, float(density), correlation


def _median_cadence_minutes(frames) -> float:
    times = [frame.midpoint_jd for frame in frames if frame.midpoint_jd is not None]
    if len(times) < 2:
        return 1.0
    differences = np.diff(np.asarray(times, dtype=float)) * 1440.0
    differences = differences[differences > 0]
    return float(np.median(differences)) if differences.size else 1.0


def _median_positive(values) -> float | None:
    accepted = [float(value) for value in values if value is not None and np.isfinite(value) and value > 0]
    return float(np.median(accepted)) if accepted else None


def _motion_limits(
    *,
    pixel_scale_arcsec: float | None,
    cadence_minutes: float,
    exposure_seconds: float,
    fwhm_px: float,
    registration_rms_px: float,
) -> tuple[float, float]:
    cadence_hours = max(cadence_minutes, 1e-3) / 60.0
    if pixel_scale_arcsec is not None:
        min_step = 8.0 / pixel_scale_arcsec * cadence_hours
        prior_max_px_hour = 600.0 / pixel_scale_arcsec
        trail_max_px_hour = 12.0 * 3600.0 / max(exposure_seconds, 1.0)
        max_step = min(prior_max_px_hour, trail_max_px_hour) * cadence_hours
    else:
        min_step = 0.25 * max(1.0, registration_rms_px)
        max_step = max(12.0, 4.0 * fwhm_px)
    min_step = float(np.clip(max(min_step, 0.18 + 0.12 * registration_rms_px), 0.18, 2.0))
    max_step = float(np.clip(max(max_step, 2.0 * fwhm_px, min_step * 4.0), 3.0, 60.0))
    return min_step, max_step
