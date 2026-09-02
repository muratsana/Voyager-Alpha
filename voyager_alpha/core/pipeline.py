from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
from PyQt6.QtCore import QThread, pyqtSignal
from astropy.wcs import WCS

from .astrometry import estimate_pixel_scale_arcsec
from .calibration import calibrate_science_frame, load_master_dark, load_master_frame
from .detection import (
    build_defect_mask,
    build_static_sky_model,
    create_residual,
    extract_residual_detections,
    robust_location_scale,
    suppress_persistent_residuals,
)
from .fits_io import read_fits_header, read_fits_image
from .known_objects import (
    KnownObjectMatcher,
    estimate_gaia_visible_limit,
    estimate_visible_limit,
    prediction_confidence,
    prediction_status,
)
from .metadata import inspect_sequence
from .models import RegistrationSolution, SequenceResult
from .plate_solver import AstapPlateSolver, merge_wcs_header
from .registration import register_frame, warp_affine
from .tracklet import link_tracklets
from .discovery_method import DOCUMENTED_DISCOVERY_METHOD


class AsteroidWorker(QThread):
    progress = pyqtSignal(int, str)
    stage_changed = pyqtSignal(str, str)
    file_progress = pyqtSignal(int, int, str, str)
    frame_done = pyqtSignal(int, dict)
    sequence_ready = pyqtSignal(list)
    known_objects_done = pyqtSignal(list)
    tracklets_done = pyqtSignal(list)
    result_ready = pyqtSignal(object)
    finished_scan = pyqtSignal()
    log = pyqtSignal(str)

    def __init__(
        self,
        fits_files,
        sigma=DOCUMENTED_DISCOVERY_METHOD.detection_sigma,
        min_pix=5,
        min_tracklet_frames=DOCUMENTED_DISCOVERY_METHOD.min_linked_frames,
        match_known_objects=True,
        master_dark_path=None,
        master_bias_path=None,
        master_flat_path=None,
        auto_plate_solve=True,
        max_sources_per_frame=DOCUMENTED_DISCOVERY_METHOD.max_residuals_per_frame,
        reference_samples=0,
        edge_margin=DOCUMENTED_DISCOVERY_METHOD.edge_margin_px,
        expected_fwhm_px=DOCUMENTED_DISCOVERY_METHOD.detection_fwhm_px,
        min_motion_px_per_frame=DOCUMENTED_DISCOVERY_METHOD.min_seed_displacement_px,
        max_step_px=35.0,
        min_median_snr=4.5,
        max_fit_rms_px=DOCUMENTED_DISCOVERY_METHOD.borderline_review_rms_px,
        strong_fit_rms_px=DOCUMENTED_DISCOVERY_METHOD.potential_discovery_rms_px,
        match_tolerance_px=DOCUMENTED_DISCOVERY_METHOD.match_tolerance_px,
        min_track_occupancy=0.5,
        max_missing_gap_frames=3,
        max_artifact_fraction=0.5,
        persistence_fraction=0.12,
        estimate_gaia_depth=True,
    ):
        super().__init__()
        self.fits_files = list(fits_files)
        self.sigma = float(sigma)
        self.min_pix = int(min_pix)
        self.min_tracklet_frames = int(min_tracklet_frames)
        self.match_known_objects = bool(match_known_objects)
        self.master_dark_path = master_dark_path or ""
        self.master_bias_path = master_bias_path or ""
        self.master_flat_path = master_flat_path or ""
        self.auto_plate_solve = bool(auto_plate_solve)
        self.max_sources_per_frame = int(max_sources_per_frame)
        self.reference_samples = max(0, int(reference_samples))
        self.edge_margin = max(4, int(edge_margin))
        self.expected_fwhm_px = float(expected_fwhm_px) if expected_fwhm_px else None
        self.min_motion_px_per_frame = max(0.05, float(min_motion_px_per_frame))
        self.max_step_px = max(self.min_motion_px_per_frame, float(max_step_px))
        self.min_median_snr = max(1.0, float(min_median_snr))
        self.max_fit_rms_px = max(0.3, float(max_fit_rms_px))
        self.strong_fit_rms_px = max(0.2, min(float(strong_fit_rms_px), self.max_fit_rms_px))
        self.match_tolerance_px = max(0.8, float(match_tolerance_px))
        self.min_track_occupancy = float(np.clip(min_track_occupancy, 0.0, 1.0))
        self.max_missing_gap_frames = max(0, int(max_missing_gap_frames))
        self.max_artifact_fraction = float(np.clip(max_artifact_fraction, 0.0, 1.0))
        self.persistence_fraction = float(np.clip(persistence_fraction, 0.05, 0.8))
        self.estimate_gaia_depth = bool(estimate_gaia_depth)
        self.is_running = True

    def run(self):
        try:
            result = self._run_pipeline()
            if result is not None and self.is_running:
                self.result_ready.emit(result)
                self.progress.emit(100, "Analiz tamamlandı")
                self.stage_changed.emit("complete", "Tamamlandı")
        except Exception as exc:
            self.log.emit(f"ERROR|Engine|{exc}")
            self.stage_changed.emit("error", str(exc))
        finally:
            self.finished_scan.emit()

    def _run_pipeline(self) -> SequenceResult | None:
        if len(self.fits_files) < 3:
            raise ValueError("Hareket analizi için en az üç FITS karesi gerekir.")

        self._stage("metadata", "Sekans doğrulanıyor", 2)
        frames = inspect_sequence(self.fits_files)
        if len(frames) < 3:
            raise ValueError("Analize uygun en az üç FITS karesi bulunamadı.")
        blocking = {
            flag
            for frame in frames
            for flag in frame.quality_flags
            if flag in {"invalid_time", "shape_mismatch"}
        }
        if blocking:
            raise ValueError(f"Engelleyici sekans hatası: {', '.join(sorted(blocking))}")
        self.fits_files = [frame.file_path for frame in frames]
        self.sequence_ready.emit(frames)
        self.log.emit(f"INFO|Metadata|{len(frames)} zaman sıralı FITS karesi doğrulandı")

        masters, dark_exposure = self._load_calibration(frames[0].shape)
        defect_mask = build_defect_mask(masters["dark"], frames[0].shape)
        reference_index = self._choose_reference(frames)
        reference_frame = frames[reference_index]
        reference_data, reference_header = self._load_calibrated(reference_frame, masters, dark_exposure)
        reference_frame.background_median, reference_frame.background_rms = robust_location_scale(reference_data)

        self._stage("solve", "Referans astrometrisi doğrulanıyor", 12)
        reference_header = self._resolve_reference_wcs(reference_frame, reference_header)
        has_reference_wcs = _has_celestial_wcs(reference_header)
        if has_reference_wcs:
            reference_frame.has_wcs = True
            reference_frame.quality_flags = [flag for flag in reference_frame.quality_flags if flag != "missing_or_invalid_WCS"]
            self.log.emit("INFO|Astrometry|Referans WCS hazır")
        else:
            self.log.emit("WARN|Astrometry|WCS üretilemedi; katalog eşleştirmesi devre dışı kalacak")

        known_predictions = []
        limiting_magnitude = None
        if has_reference_wcs and self.match_known_objects:
            self._stage("known", "Generate: SkyBoT tahminleri ve yerel piksel kanıtı üretiliyor", 16)
            matcher = KnownObjectMatcher(tolerance_arcsec=20.0)
            known_predictions = matcher.predictions_for_frame(
                reference_frame,
                reference_header,
                image=reference_data,
            )
            limiting_magnitude = None
            if self.estimate_gaia_depth:
                self.log.emit("INFO|Estimate|Gaia G yıldızlarıyla görünürlük limiti ölçülüyor")
                limiting_magnitude = estimate_gaia_visible_limit(reference_data, reference_header)
            if limiting_magnitude is None:
                limiting_magnitude = estimate_visible_limit(known_predictions)
            if limiting_magnitude is not None:
                for prediction in known_predictions:
                    prediction.confidence = prediction_confidence(prediction, limiting_magnitude)
                    prediction.status = prediction_status(prediction)
            self.known_objects_done.emit(known_predictions)
            self.log.emit(
                f"INFO|Generate|{len(known_predictions)} alan tahmini; "
                f"{sum(1 for item in known_predictions if item.visible)} yerel piksel eşleşmesi"
            )
        elif self.match_known_objects:
            matcher = None
            self.log.emit("WARN|Generate|WCS olmadığı için bilinen cisim tahmini üretilemedi")
        else:
            matcher = None

        self._stage("reference", "Align: tüm kareler ortak piksel ızgarasına getiriliyor", 18)
        registration_solutions: dict[int, RegistrationSolution] = {
            reference_index: RegistrationSolution.identity()
        }
        detections_by_frame = {}
        reliable_registration_positions: set[int] = set()
        with tempfile.TemporaryDirectory(prefix="astrohub-discover-") as cache_dir:
            cache_path = Path(cache_dir) / "aligned.float32"
            aligned_stack = np.memmap(
                cache_path,
                dtype=np.float32,
                mode="w+",
                shape=(len(frames), *reference_data.shape),
            )
            for position, frame in enumerate(frames):
                if not self.is_running:
                    return None
                frame_index = frame.index
                self.file_progress.emit(position + 1, len(frames), Path(frame.file_path).name, "Align")
                if frame_index == reference_index:
                    aligned = reference_data
                    solution = registration_solutions[reference_index]
                    reliable_registration_positions.add(position)
                else:
                    data, _header = self._load_calibrated(frame, masters, dark_exposure)
                    aligned, solution = register_frame(reference_data, data)
                    registration_solutions[frame_index] = solution
                    self._store_registration(frame, solution)
                    if "rot180" in solution.method:
                        self.log.emit(
                            f"INFO|Align|{Path(frame.file_path).name}: meridian flip algılandı; "
                            "180 derece yön düzeltmesi ve yıldız affine hizası uygulandı"
                        )
                    if solution.method.startswith("star") and solution.rms_px <= 2.0:
                        reliable_registration_positions.add(position)
                    else:
                        if "weak_star_registration" not in frame.quality_flags:
                            frame.quality_flags.append("weak_star_registration")
                        self.log.emit(
                            f"WARN|Align|{Path(frame.file_path).name}: güvenilir yıldız hizası yok "
                            f"({solution.method}, {solution.matched_stars} eşleşme); aday üretiminden çıkarıldı"
                        )
                aligned_stack[position] = np.asarray(aligned, dtype=np.float32)
                self.progress.emit(18 + int(18 * (position + 1) / len(frames)), f"Align {position + 1}/{len(frames)}")
            aligned_stack.flush()
            if len(reliable_registration_positions) < 3:
                raise ValueError(
                    "En az üç karede güvenilir yıldız hizalaması kurulamadı. "
                    "Tespit çalıştırılmadı; eşik, pozlama ve yıldız görünürlüğünü kontrol edin."
                )
            static_model = _temporal_median_from_memmap(
                aligned_stack,
                frame_indices=sorted(reliable_registration_positions),
            )
            self.log.emit("INFO|Reference|Tüm hizalı karelerden exact temporal-median statik gökyüzü modeli üretildi")

            self._stage("detect", "Discover: hibrit nokta/iz residual taraması", 38)
            for position, frame in enumerate(frames):
                if not self.is_running:
                    return None
                frame_index = frame.index
                name = Path(frame.file_path).name
                self.file_progress.emit(position + 1, len(frames), name, "Discover")
                aligned = np.asarray(aligned_stack[position], dtype=np.float32)
                solution = registration_solutions[frame_index]
                residual = create_residual(aligned, static_model)
                if position in reliable_registration_positions:
                    detections, residual_rms = extract_residual_detections(
                        residual,
                        frame_index,
                        header=reference_header if has_reference_wcs else None,
                        threshold_sigma=self.sigma,
                        min_pixels=self.min_pix,
                        max_sources=self.max_sources_per_frame,
                        edge_margin=self.edge_margin,
                        expected_fwhm_px=self.expected_fwhm_px,
                        defect_mask=defect_mask,
                        detector_mode=DOCUMENTED_DISCOVERY_METHOD.detector_mode,
                        streak_min_area=DOCUMENTED_DISCOVERY_METHOD.streak_min_area_px,
                        streak_min_elongation=DOCUMENTED_DISCOVERY_METHOD.streak_min_elongation,
                    )
                else:
                    detections = []
                    _median, residual_rms = robust_location_scale(residual)
                for detection in detections:
                    detection.sensor_x, detection.sensor_y = _aligned_to_sensor_xy(
                        detection.x,
                        detection.y,
                        solution,
                    )
                frame.background_median, frame.background_rms = robust_location_scale(aligned)
                detections_by_frame[frame_index] = detections
                preview = _preview_image(residual)
                self.frame_done.emit(
                    frame_index,
                    {
                        "file": frame.file_path,
                        "diff_data": preview,
                        "residual_preview": preview,
                        "candidates": [detection.to_overlay() for detection in detections],
                        "detections": detections,
                        "registration": solution,
                        "residual_rms": residual_rms,
                    },
                )
                self.progress.emit(38 + int(38 * (position + 1) / len(frames)), f"Discover {position + 1}/{len(frames)}")
                self.log.emit(f"INFO|Detect|{name}: {len(detections)} residual kaynak")
            del aligned
            del aligned_stack

        raw_detection_count = sum(len(items) for items in detections_by_frame.values())
        detections_by_frame = suppress_persistent_residuals(
            detections_by_frame,
            persistence_fraction=self.persistence_fraction,
        )
        filtered_detection_count = sum(len(items) for items in detections_by_frame.values())
        rejected_artifacts = raw_detection_count - filtered_detection_count
        self.log.emit(
            f"INFO|Artifact filter|{rejected_artifacts} sıcak piksel, sensör-sabit kusur veya kalıcı residual elendi; "
            f"{filtered_detection_count} tek-kare residual bağlantı aşamasına geçti"
        )
        self._stage("link", "Zaman tabanlı tracklet modeli kuruluyor", 78)
        pixel_scale = estimate_pixel_scale_arcsec(reference_header) if has_reference_wcs else None
        tracklets = link_tracklets(
            detections_by_frame,
            frames,
            min_frames=self.min_tracklet_frames,
            min_median_snr=max(self.min_median_snr, self.sigma * 0.9),
            min_motion_px_per_frame=self.min_motion_px_per_frame,
            max_step_px=self.max_step_px,
            max_fit_rms_px=self.max_fit_rms_px,
            strong_fit_rms_px=self.strong_fit_rms_px,
            match_tolerance_px=self.match_tolerance_px,
            min_track_occupancy=self.min_track_occupancy,
            max_missing_gap_frames=self.max_missing_gap_frames,
            max_artifact_fraction=self.max_artifact_fraction,
            pixel_scale_arcsec=pixel_scale,
        )

        if matcher is not None:
            self._stage("known", "Recover: trackletler katalog tahminleriyle karşılaştırılıyor", 88)
            for tracklet in tracklets:
                representative = tracklet.detections[len(tracklet.detections) // 2]
                frame = frames[representative.frame_index]
                match = matcher.match_tracklet(tracklet, frame, reference_header)
                if match:
                    tracklet.known_match = match
                    tracklet.classification = "known_object"
            self.log.emit(
                f"INFO|Catalog|{len(known_predictions)} bilinen cisim tahmini, "
                f"{sum(1 for item in known_predictions if item.visible)} yerel görüntü eşleşmesi"
            )

        self.tracklets_done.emit(tracklets)
        unknown_count = sum(1 for item in tracklets if item.classification in {"unknown_candidate", "review_candidate"})
        review_count = sum(1 for item in tracklets if item.classification == "review_candidate")
        level = "WARN" if unknown_count else "INFO"
        self.log.emit(
            f"{level}|Classify|{len(tracklets)} hareket zinciri; {unknown_count} bilinmeyen aday, "
            f"{review_count} sınırda inceleme. "
            "Tüm adaylar insan doğrulaması gerektirir"
        )
        return SequenceResult(
            frames=frames,
            reference_index=reference_index,
            reference_header=reference_header,
            known_objects=known_predictions,
            tracklets=tracklets,
            limiting_magnitude=limiting_magnitude,
            registration_solutions=registration_solutions,
            static_model=static_model,
            reference_data=reference_data,
            calibration_paths={
                key: value
                for key, value in {
                    "dark": self.master_dark_path,
                    "bias": self.master_bias_path,
                    "flat": self.master_flat_path,
                }.items()
                if value
            },
            detections_by_frame=detections_by_frame,
        )

    def _load_calibration(self, expected_shape):
        masters = {"bias": None, "dark": None, "flat": None}
        dark_exposure = None
        if self.master_bias_path:
            masters["bias"] = load_master_frame(self.master_bias_path)
        if self.master_dark_path:
            masters["dark"], dark_exposure = load_master_dark(self.master_dark_path)
        if self.master_flat_path:
            masters["flat"] = load_master_frame(self.master_flat_path)
        for name, master in masters.items():
            if master is not None and master.shape != expected_shape:
                raise ValueError(f"Master {name} boyutu sekansla uyuşmuyor: {master.shape} != {expected_shape}")
        active = [name for name, value in masters.items() if value is not None]
        self.log.emit(
            f"INFO|Calibration|Uygulanan master kareler: {', '.join(active) if active else 'yok'}"
        )
        return masters, dark_exposure

    def _load_calibrated(self, frame, masters, dark_exposure):
        data, header = read_fits_image(frame.file_path, header=True)
        calibrated = calibrate_science_frame(
            data,
            master_bias=masters["bias"],
            master_dark=masters["dark"],
            master_flat=masters["flat"],
            science_exposure=frame.exposure_seconds,
            dark_exposure=dark_exposure,
        )
        frame.calibration_state = "+".join(name for name, value in masters.items() if value is not None) or "raw"
        return calibrated, header

    def _choose_reference(self, frames):
        center = (len(frames) - 1) / 2.0
        return min(
            range(len(frames)),
            key=lambda index: (0 if frames[index].has_wcs else 1, abs(index - center), len(frames[index].quality_flags)),
        )

    def _resolve_reference_wcs(self, frame, header):
        if _has_celestial_wcs(header):
            return header
        if not self.auto_plate_solve:
            return header
        solver = AstapPlateSolver()
        if not solver.is_available():
            self.log.emit("WARN|Astrometry|ASTAP bulunamadı")
            return header
        self.log.emit(f"INFO|Astrometry|ASTAP referans çözümü başladı: {Path(frame.file_path).name}")
        result = solver.solve(frame.file_path, fov_deg=0.0)
        if result.success and result.header is not None:
            return merge_wcs_header(header, result.header)
        self.log.emit(f"WARN|Astrometry|ASTAP çözümü başarısız: {result.message[:300]}")
        return header

    def _store_registration(self, frame, solution):
        frame.registration_dx = float(solution.offset_xy[0])
        frame.registration_dy = float(solution.offset_xy[1])
        frame.registration_peak = float(solution.phase_peak)
        frame.registration_rms_px = float(solution.rms_px)
        frame.matched_stars = int(solution.matched_stars)
        if solution.rms_px > 2.0:
            frame.quality_flags.append("high_registration_rms")
        if not solution.method.startswith("star"):
            frame.quality_flags.append("weak_star_affine_registration")

    def _stage(self, key, label, progress):
        self.stage_changed.emit(key, label)
        self.progress.emit(int(progress), label)
        self.log.emit(f"INFO|Workflow|{label}")

    def stop(self):
        self.is_running = False


def _aligned_to_sensor_xy(x: float, y: float, solution: RegistrationSolution) -> tuple[float, float]:
    aligned = np.asarray([float(x), float(y)], dtype=np.float64)
    try:
        sensor = np.linalg.solve(solution.matrix_xy, aligned - solution.offset_xy)
        if np.all(np.isfinite(sensor)):
            return float(sensor[0]), float(sensor[1])
    except (ValueError, np.linalg.LinAlgError):
        pass
    return float(x), float(y)


def materialize_sequence_frame(
    result: SequenceResult,
    frame_index: int,
    *,
    residual: bool = False,
) -> np.ndarray:
    if frame_index < 0 or frame_index >= len(result.frames):
        raise IndexError(frame_index)
    frame = result.frames[frame_index]
    masters = {"bias": None, "dark": None, "flat": None}
    dark_exposure = None
    if result.calibration_paths.get("bias"):
        masters["bias"] = load_master_frame(result.calibration_paths["bias"])
    if result.calibration_paths.get("dark"):
        masters["dark"], dark_exposure = load_master_dark(result.calibration_paths["dark"])
    if result.calibration_paths.get("flat"):
        masters["flat"] = load_master_frame(result.calibration_paths["flat"])
    data = calibrate_science_frame(
        read_fits_image(frame.file_path),
        master_bias=masters["bias"],
        master_dark=masters["dark"],
        master_flat=masters["flat"],
        science_exposure=frame.exposure_seconds,
        dark_exposure=dark_exposure,
    )
    if frame_index == result.reference_index:
        aligned = data
    else:
        solution = result.registration_solutions[frame_index]
        aligned = warp_affine(
            data,
            solution,
            result.reference_data.shape,
            fill_value=float(np.nanmedian(result.reference_data)),
        )
    if residual:
        return create_residual(aligned, result.static_model)
    return aligned


def _preview_image(image: np.ndarray, max_dim: int = 1400) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    factor = max(1, int(np.ceil(max(data.shape) / max_dim)))
    return data[::factor, ::factor].astype(np.float32, copy=False)


def _temporal_median_from_memmap(
    aligned_stack: np.memmap,
    *,
    target_working_mb: int = 192,
    frame_indices: list[int] | None = None,
) -> np.ndarray:
    """Compute an exact temporal median in row bands without loading the cube into RAM."""

    frame_count, height, width = aligned_stack.shape
    selected = list(range(frame_count)) if frame_indices is None else [int(index) for index in frame_indices]
    if not selected:
        raise ValueError("Temporal median requires at least one aligned frame.")
    bytes_per_row = max(1, len(selected) * width * np.dtype(np.float32).itemsize)
    rows_per_band = max(1, min(height, target_working_mb * 1024 * 1024 // bytes_per_row))
    model = np.empty((height, width), dtype=np.float32)
    for row_start in range(0, height, rows_per_band):
        row_end = min(height, row_start + rows_per_band)
        band = np.asarray(aligned_stack[selected, row_start:row_end, :], dtype=np.float32)
        model[row_start:row_end] = np.nanmedian(band, axis=0).astype(np.float32)
    return model


def _has_celestial_wcs(header) -> bool:
    try:
        return bool(WCS(header).has_celestial)
    except Exception:
        return False
