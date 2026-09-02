from __future__ import annotations

import numpy as np
import sep
from astropy.wcs import WCS
from scipy import ndimage

from .models import Detection
from .discovery_method import DOCUMENTED_DISCOVERY_METHOD


def robust_location_scale(data: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(data, dtype=np.float32)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    median = float(np.median(finite))
    mad = float(np.median(np.abs(finite - median)))
    scale = max(1.4826 * mad, float(np.std(finite)) * 0.08, 1e-6)
    return median, scale


def build_static_sky_model(aligned_samples: list[np.ndarray]) -> np.ndarray:
    if len(aligned_samples) < 3:
        raise ValueError("Static sky model requires at least three aligned frames.")
    shape = aligned_samples[0].shape
    normalized = []
    for frame in aligned_samples:
        if frame.shape != shape:
            raise ValueError("Aligned frame shapes do not match.")
        normalized.append(np.asarray(frame, dtype=np.float32))
    cube = np.stack(normalized, axis=0)
    model = np.nanmedian(cube, axis=0).astype(np.float32)
    del cube
    return model


def match_background_and_scale(image: np.ndarray, model: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    reference = np.asarray(model, dtype=np.float32)
    step = max(1, int(np.ceil(max(data.shape) / 700)))
    x = reference[::step, ::step].ravel().astype(np.float64)
    y = data[::step, ::step].ravel().astype(np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < 100:
        return data.copy()
    x_median, x_scale = robust_location_scale(x)
    y_median, y_scale = robust_location_scale(y)
    keep = (
        (np.abs(x - x_median) < 4.0 * x_scale)
        & (np.abs(y - y_median) < 4.0 * y_scale)
    )
    if int(np.count_nonzero(keep)) >= 100:
        design = np.column_stack((y[keep], np.ones(int(np.count_nonzero(keep)))))
        slope, intercept = np.linalg.lstsq(design, x[keep], rcond=None)[0]
        if np.isfinite(slope) and 0.25 <= slope <= 4.0:
            return (data * float(slope) + float(intercept)).astype(np.float32)
    scale = x_scale / max(y_scale, 1e-6)
    return ((data - y_median) * scale + x_median).astype(np.float32)


def create_residual(image: np.ndarray, static_model: np.ndarray) -> np.ndarray:
    matched = match_background_and_scale(image, static_model)
    return (matched - np.asarray(static_model, dtype=np.float32)).astype(np.float32)


def extract_residual_detections(
    residual: np.ndarray,
    frame_index: int,
    *,
    header=None,
    threshold_sigma: float = 5.0,
    min_pixels: int = 5,
    max_sources: int = DOCUMENTED_DISCOVERY_METHOD.max_residuals_per_frame,
    edge_margin: int = DOCUMENTED_DISCOVERY_METHOD.edge_margin_px,
    expected_fwhm_px: float | None = None,
    defect_mask: np.ndarray | None = None,
    max_elongation: float = 6.0,
    detector_mode: str = DOCUMENTED_DISCOVERY_METHOD.detector_mode,
    streak_min_area: int = DOCUMENTED_DISCOVERY_METHOD.streak_min_area_px,
    streak_min_elongation: float = DOCUMENTED_DISCOVERY_METHOD.streak_min_elongation,
) -> tuple[list[Detection], float]:
    if detector_mode not in {"point", "streak", "hybrid"}:
        raise ValueError(f"Unsupported residual detector mode: {detector_mode}")
    data = _native_float32(residual)
    background = sep.Background(data)
    centered = data - background
    _median, robust_rms = robust_location_scale(centered)
    effective_rms = max(float(background.globalrms), robust_rms, 1e-6)
    objects = sep.extract(
        centered,
        float(threshold_sigma),
        err=effective_rms,
        minarea=int(min_pixels),
        deblend_nthresh=32,
        deblend_cont=0.01,
    )
    if len(objects) == 0:
        return [], effective_rms
    height, width = data.shape
    detections = []
    if defect_mask is not None and np.asarray(defect_mask).shape != data.shape:
        raise ValueError("Defect mask shape does not match residual frame.")
    wcs = None
    if header is not None:
        try:
            candidate_wcs = WCS(header)
            wcs = candidate_wcs if candidate_wcs.has_celestial else None
        except Exception:
            wcs = None

    for object_index in range(len(objects)):
        obj = objects[object_index]
        x, y = float(obj["x"]), float(obj["y"])
        if x < edge_margin or y < edge_margin or x >= width - edge_margin or y >= height - edge_margin:
            continue
        major = max(float(obj["a"]), 1e-3)
        minor = max(float(obj["b"]), 1e-3)
        ratio = major / minor
        measured_fwhm = float(2.355 * np.sqrt(major * minor))
        minor_fwhm = float(2.355 * minor)
        area = int(obj["npix"]) if "npix" in obj.dtype.names else int(min_pixels)
        flux = float(obj["flux"])
        peak = float(obj["peak"]) if "peak" in obj.dtype.names else flux / max(area, 1)
        local_rms = max(effective_rms, _local_background_rms(centered, x, y, measured_fwhm))
        snr = flux / max(local_rms * np.sqrt(max(area, 1)), 1e-6)
        support_pixels, effective_area = _positive_source_support(
            centered,
            x,
            y,
            local_rms,
            peak,
            expected_fwhm_px,
        )
        flags = []
        if defect_mask is not None and bool(defect_mask[int(round(y)), int(round(x))]):
            continue
        if ratio > max_elongation:
            continue
        is_streak = area >= int(streak_min_area) and ratio >= float(streak_min_elongation)
        if detector_mode == "point" and is_streak:
            continue
        if detector_mode == "streak" and not is_streak:
            continue
        flags.append("streak_source" if is_streak else "point_source")
        if ratio > 3.2:
            flags.append("elongated_or_trail")
        minimum_fwhm = max(1.25, float(expected_fwhm_px) * 0.52) if expected_fwhm_px is not None else 1.25
        minimum_minor_fwhm = max(1.05, float(expected_fwhm_px) * 0.42) if expected_fwhm_px is not None else 1.05
        if measured_fwhm < minimum_fwhm or minor_fwhm < minimum_minor_fwhm:
            flags.append("hot_pixel_or_cosmic_ray")
            continue
        if area < max(3, int(min_pixels)):
            flags.append("hot_pixel_or_cosmic_ray")
            continue
        if support_pixels < max(5, int(min_pixels)) or effective_area < 1.8:
            flags.append("hot_pixel_or_cosmic_ray")
            continue
        if _looks_saturated(centered, x, y, peak):
            flags.append("near_saturation")
            saturation_fwhm_limit = max(2.2, float(expected_fwhm_px) * 0.8) if expected_fwhm_px else 2.2
            if measured_fwhm < saturation_fwhm_limit:
                continue
        if _has_negative_dipole(centered, x, y, local_rms, threshold_sigma):
            flags.append("star_subtraction_dipole")
        if peak / max(local_rms, 1e-6) < float(threshold_sigma) * 0.75:
            continue
        if "star_subtraction_dipole" in flags and snr < max(9.0, float(threshold_sigma) * 1.8):
            continue
        ra = dec = None
        if wcs is not None:
            try:
                ra_value, dec_value = wcs.all_pix2world(x, y, 0)
                ra, dec = float(ra_value), float(dec_value)
            except Exception:
                pass
        detections.append(
            Detection(
                frame_index=int(frame_index),
                x=x,
                y=y,
                ra=ra,
                dec=dec,
                flux=flux,
                snr=float(snr),
                fwhm_px=measured_fwhm,
                eccentricity=float(np.sqrt(max(0.0, 1.0 - (minor / major) ** 2))),
                area_px=area,
                peak=peak,
                local_rms=local_rms,
                flags=flags,
            )
        )
    detections.sort(key=_detection_rank, reverse=True)
    return detections[: max(1, int(max_sources))], effective_rms


def suppress_persistent_residuals(
    detections_by_frame: dict[int, list[Detection]],
    *,
    cell_size_px: float = 2.5,
    persistence_fraction: float = 0.12,
) -> dict[int, list[Detection]]:
    frame_count = max(1, len(detections_by_frame))
    minimum_hits = max(3, int(np.ceil(frame_count * persistence_fraction)))
    cell_frames: dict[tuple[int, int], set[int]] = {}
    cell_positions: dict[tuple[int, int], list[tuple[int, float, float]]] = {}
    sensor_cell_frames: dict[tuple[int, int], set[int]] = {}
    sensor_cell_positions: dict[tuple[int, int], list[tuple[float, float]]] = {}
    for frame_index, detections in detections_by_frame.items():
        for detection in detections:
            cell = (int(round(detection.x / cell_size_px)), int(round(detection.y / cell_size_px)))
            cell_frames.setdefault(cell, set()).add(frame_index)
            cell_positions.setdefault(cell, []).append((frame_index, detection.x, detection.y))
            if detection.sensor_x is not None and detection.sensor_y is not None:
                sensor_cell = (
                    int(round(detection.sensor_x / cell_size_px)),
                    int(round(detection.sensor_y / cell_size_px)),
                )
                sensor_cell_frames.setdefault(sensor_cell, set()).add(frame_index)
                sensor_cell_positions.setdefault(sensor_cell, []).append((detection.sensor_x, detection.sensor_y))
    persistent = set()
    for cell, seen in cell_frames.items():
        if len(seen) < minimum_hits:
            continue
        samples = np.asarray(cell_positions[cell], dtype=float)
        frame_samples = samples[:, 0]
        points = samples[:, 1:]
        center = np.median(points, axis=0)
        radial_spread = np.percentile(np.hypot(points[:, 0] - center[0], points[:, 1] - center[1]), 90)
        frame_span = float(np.ptp(frame_samples))
        fitted_displacement = 0.0
        if frame_span > 0.0:
            relative_frames = frame_samples - frame_samples.min()
            velocity_x = float(np.polyfit(relative_frames, points[:, 0], 1)[0])
            velocity_y = float(np.polyfit(relative_frames, points[:, 1], 1)[0])
            fitted_displacement = float(np.hypot(velocity_x, velocity_y) * frame_span)
        if radial_spread <= 0.70 and fitted_displacement <= 0.55:
            persistent.add(cell)
    persistent_sensor = set()
    for cell, seen in sensor_cell_frames.items():
        if len(seen) < minimum_hits:
            continue
        points = np.asarray(sensor_cell_positions[cell], dtype=float)
        center = np.median(points, axis=0)
        radial_spread = np.percentile(np.hypot(points[:, 0] - center[0], points[:, 1] - center[1]), 90)
        if radial_spread <= 0.70:
            persistent_sensor.add(cell)
    filtered = {}
    for frame_index, detections in detections_by_frame.items():
        accepted = []
        for detection in detections:
            cell = (int(round(detection.x / cell_size_px)), int(round(detection.y / cell_size_px)))
            if cell in persistent:
                detection.flags.append("persistent_star_residual")
                continue
            if detection.sensor_x is not None and detection.sensor_y is not None:
                sensor_cell = (
                    int(round(detection.sensor_x / cell_size_px)),
                    int(round(detection.sensor_y / cell_size_px)),
                )
                if sensor_cell in persistent_sensor:
                    detection.flags.append("persistent_sensor_defect")
                    continue
            if "star_subtraction_dipole" in detection.flags and detection.snr < 9.0:
                continue
            accepted.append(detection)
        filtered[frame_index] = accepted
    return filtered


def build_defect_mask(
    master_dark: np.ndarray | None,
    shape: tuple[int, int],
) -> np.ndarray:
    if master_dark is None:
        return np.zeros(shape, dtype=bool)
    dark = np.asarray(master_dark, dtype=np.float32)
    if dark.shape != shape:
        raise ValueError("Master dark defect mask shape mismatch.")
    median, scale = robust_location_scale(dark)
    hot = dark > median + 8.0 * scale
    return ndimage.binary_dilation(hot, iterations=1)


def _has_negative_dipole(
    residual: np.ndarray,
    x: float,
    y: float,
    rms: float,
    threshold_sigma: float,
) -> bool:
    radius = 6
    x0, x1 = max(0, int(x) - radius), min(residual.shape[1], int(x) + radius + 1)
    y0, y1 = max(0, int(y) - radius), min(residual.shape[0], int(y) + radius + 1)
    patch = residual[y0:y1, x0:x1]
    if patch.size == 0:
        return False
    negative = abs(min(float(np.nanmin(patch)), 0.0))
    positive = max(float(np.nanmax(patch)), 0.0)
    return negative >= threshold_sigma * rms * 0.75 and negative >= positive * 0.35


def _local_background_rms(residual: np.ndarray, x: float, y: float, fwhm_px: float) -> float:
    radius = int(np.clip(round(max(8.0, 3.0 * fwhm_px)), 8, 24))
    x0, x1 = max(0, int(x) - radius), min(residual.shape[1], int(x) + radius + 1)
    y0, y1 = max(0, int(y) - radius), min(residual.shape[0], int(y) + radius + 1)
    patch = residual[y0:y1, x0:x1]
    if patch.size < 25:
        return 0.0
    yy, xx = np.indices(patch.shape)
    center_x = x - x0
    center_y = y - y0
    exclusion = np.hypot(xx - center_x, yy - center_y) <= max(2.5, 1.5 * fwhm_px)
    background_pixels = patch[~exclusion]
    _median, scale = robust_location_scale(background_pixels)
    return float(scale)


def _positive_source_support(
    residual: np.ndarray,
    x: float,
    y: float,
    rms: float,
    peak: float,
    expected_fwhm_px: float | None,
) -> tuple[int, float]:
    radius = int(np.clip(np.ceil(1.5 * (expected_fwhm_px or 2.5)), 2, 7))
    center_x = int(round(x))
    center_y = int(round(y))
    x0, x1 = max(0, center_x - radius), min(residual.shape[1], center_x + radius + 1)
    y0, y1 = max(0, center_y - radius), min(residual.shape[0], center_y + radius + 1)
    patch = np.asarray(residual[y0:y1, x0:x1], dtype=np.float32)
    if patch.size == 0 or not np.isfinite(peak) or peak <= 0:
        return 0, 0.0
    local_x = int(np.clip(center_x - x0, 0, patch.shape[1] - 1))
    local_y = int(np.clip(center_y - y0, 0, patch.shape[0] - 1))
    search = patch[
        max(0, local_y - 1) : min(patch.shape[0], local_y + 2),
        max(0, local_x - 1) : min(patch.shape[1], local_x + 2),
    ]
    if search.size:
        peak_offset = np.unravel_index(int(np.nanargmax(search)), search.shape)
        local_y = max(0, local_y - 1) + int(peak_offset[0])
        local_x = max(0, local_x - 1) + int(peak_offset[1])
    threshold = max(1.5 * max(float(rms), 1e-6), 0.08 * float(peak))
    mask = np.isfinite(patch) & (patch > threshold)
    labels, _count = ndimage.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    label = int(labels[local_y, local_x])
    if label == 0:
        return 0, 0.0
    component = labels == label
    support = int(np.count_nonzero(component))
    positive_flux = float(np.sum(np.clip(patch[component], 0.0, None)))
    return support, positive_flux / max(float(peak), 1e-6)


def _looks_saturated(residual: np.ndarray, x: float, y: float, peak: float) -> bool:
    if not np.isfinite(peak) or peak <= 0:
        return False
    radius = 3
    x0, x1 = max(0, int(x) - radius), min(residual.shape[1], int(x) + radius + 1)
    y0, y1 = max(0, int(y) - radius), min(residual.shape[0], int(y) + radius + 1)
    patch = residual[y0:y1, x0:x1]
    return int(np.count_nonzero(patch >= peak * 0.985)) >= 4


def _detection_rank(detection: Detection) -> float:
    penalty = 0.0
    if "star_subtraction_dipole" in detection.flags:
        penalty += 4.0
    if "near_saturation" in detection.flags:
        penalty += 2.0
    return float(detection.snr - penalty)


def _native_float32(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    if not data.dtype.isnative:
        data = data.byteswap().view(data.dtype.newbyteorder("="))
    return np.ascontiguousarray(data)
