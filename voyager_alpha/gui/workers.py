from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QThread, Qt, pyqtSignal
from PyQt6.QtGui import QImage
from scipy import ndimage

from ..core.fits_io import read_fits_image
from ..core.detection import robust_location_scale
from ..core.detection_optimizer import optimize_detection_settings
from ..core.models import SyntheticTrackResult
from ..core.pipeline import materialize_sequence_frame
from ..core.plate_solver import AstapPlateSolver
from ..core.registration import register_frame
from ..core.synthetic_tracking import _stack_peak_metrics, combine_centered_cutouts
from .viewer import stretch_image


class PreviewCacheWorker(QThread):
    preview_ready = pyqtSignal(int, object, object, object)
    alignment_ready = pyqtSignal(int, str, float, int)
    failed = pyqtSignal(str)
    complete = pyqtSignal()

    def __init__(self, file_paths, *, max_dim=1280, sequence_result=None, parent=None):
        super().__init__(parent)
        self.file_paths = list(file_paths)
        self.max_dim = int(max_dim)
        self.sequence_result = sequence_result
        self._stopping = False

    def stop(self):
        self._stopping = True

    def run(self):
        try:
            if not self.file_paths:
                return
            reference_preview = None
            for index, file_path in enumerate(self.file_paths):
                if self._stopping:
                    return
                try:
                    if self.sequence_result is not None:
                        data = np.asarray(
                            materialize_sequence_frame(self.sequence_result, index, residual=False),
                            dtype=np.float32,
                        )
                        preview = _downsample_preview(data, self.max_dim)
                        solution = self.sequence_result.registration_solutions.get(index)
                    else:
                        data = np.asarray(read_fits_image(file_path), dtype=np.float32)
                        preview = _downsample_preview(data, self.max_dim)
                        if reference_preview is None:
                            reference_preview = preview
                            solution = None
                        else:
                            preview, solution = register_frame(
                                reference_preview,
                                preview,
                                max_shift_px=min(180.0, 0.35 * max(reference_preview.shape)),
                            )
                    if solution is None:
                        self.alignment_ready.emit(index, "reference", 0.0, 0)
                    else:
                        self.alignment_ready.emit(
                            index,
                            solution.method,
                            float(solution.rms_px),
                            int(solution.matched_stars),
                        )
                    display, _ = stretch_image(preview, mode="Auto STF", sequence_levels=None)
                    display_u8 = np.ascontiguousarray(np.clip(display * 255.0, 0, 255), dtype=np.uint8)
                    thumbnail = _thumbnail_image(display_u8)
                    self.preview_ready.emit(index, display_u8, data.shape, thumbnail)
                except Exception as exc:
                    self.failed.emit(f"{file_path}: {exc}")
        except Exception as exc:
            self.failed.emit(str(exc))
        finally:
            self.complete.emit()


class DetectionOptimizationWorker(QThread):
    progress = pyqtSignal(int, str)
    result_ready = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, file_paths, master_paths, parent=None):
        super().__init__(parent)
        self.file_paths = list(file_paths)
        self.master_paths = dict(master_paths)
        self._stopping = False

    def stop(self):
        self._stopping = True

    def run(self):
        try:
            result = optimize_detection_settings(
                self.file_paths,
                master_paths=self.master_paths,
                progress_callback=lambda value, message: self.progress.emit(value, message),
                stop_callback=lambda: self._stopping,
            )
            if not self._stopping:
                self.result_ready.emit(result)
        except Exception as exc:
            if not self._stopping:
                self.failed.emit(str(exc))


class SyntheticTrackWorker(QThread):
    progress = pyqtSignal(int, int)
    result_ready = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, sequence_result, tracklet, *, radius_px=28, parent=None):
        super().__init__(parent)
        self.sequence_result = sequence_result
        self.tracklet = tracklet
        self.radius_px = int(radius_px)
        self._stopping = False

    def stop(self):
        self._stopping = True

    def run(self):
        try:
            detections = sorted(self.tracklet.detections, key=lambda item: item.frame_index)
            anchor = detections[len(detections) // 2]
            anchor_time = _time_minutes(self.sequence_result.frames, anchor.frame_index)
            axis = np.arange(-self.radius_px, self.radius_px + 1, dtype=np.float64)
            yy, xx = np.meshgrid(axis, axis, indexing="ij")
            cutouts = []
            skipped = 0
            total = len(self.sequence_result.frames)
            for index, frame_record in enumerate(self.sequence_result.frames):
                if self._stopping:
                    return
                aligned = materialize_sequence_frame(self.sequence_result, index, residual=False)
                dt = _time_minutes(self.sequence_result.frames, index) - anchor_time
                predicted_x = anchor.x + self.tracklet.velocity_x_px_min * dt
                predicted_y = anchor.y + self.tracklet.velocity_y_px_min * dt
                if (
                    predicted_x - self.radius_px < 0
                    or predicted_y - self.radius_px < 0
                    or predicted_x + self.radius_px >= aligned.shape[1]
                    or predicted_y + self.radius_px >= aligned.shape[0]
                ):
                    skipped += 1
                else:
                    sampled = ndimage.map_coordinates(
                        aligned,
                        [predicted_y + yy, predicted_x + xx],
                        order=1,
                        mode="nearest",
                        prefilter=False,
                    )
                    median, _rms = robust_location_scale(sampled)
                    cutouts.append((sampled - median).astype(np.float32))
                self.progress.emit(index + 1, total)
            stack = combine_centered_cutouts(cutouts)
            snr, offset = _stack_peak_metrics(stack)
            self.result_ready.emit(
                SyntheticTrackResult(
                    image=stack,
                    used_frames=len(cutouts),
                    skipped_frames=skipped,
                    velocity_x_px_min=self.tracklet.velocity_x_px_min,
                    velocity_y_px_min=self.tracklet.velocity_y_px_min,
                    snr=snr,
                    peak_offset_px=offset,
                    center_x=anchor.x,
                    center_y=anchor.y,
                )
            )
        except Exception as exc:
            self.failed.emit(str(exc))


class PlateSolveWorker(QThread):
    result_ready = pyqtSignal(object)

    def __init__(self, file_path: str, parent=None):
        super().__init__(parent)
        self.file_path = file_path

    def run(self):
        self.result_ready.emit(AstapPlateSolver().solve(self.file_path, fov_deg=0.0))


def _thumbnail_image(display_u8: np.ndarray) -> QImage:
    height, width = display_u8.shape
    image = QImage(
        display_u8.data,
        width,
        height,
        width,
        QImage.Format.Format_Grayscale8,
    ).copy()
    return image.scaled(
        126,
        72,
        aspectRatioMode=Qt.AspectRatioMode.KeepAspectRatioByExpanding,
        transformMode=Qt.TransformationMode.SmoothTransformation,
    )


def _downsample_preview(data: np.ndarray, max_dim: int) -> np.ndarray:
    factor = max(1.0, max(data.shape) / max(1, int(max_dim)))
    if factor <= 1.0:
        return np.asarray(data, dtype=np.float32)
    return ndimage.zoom(
        data,
        zoom=(1.0 / factor, 1.0 / factor),
        order=1,
        prefilter=False,
    ).astype(np.float32)


def _time_minutes(frames, index: int) -> float:
    valid = [frame.midpoint_jd for frame in frames if frame.midpoint_jd is not None]
    if not valid or frames[index].midpoint_jd is None:
        return float(index)
    return (float(frames[index].midpoint_jd) - float(valid[0])) * 1440.0
