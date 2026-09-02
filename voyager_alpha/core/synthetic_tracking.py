from __future__ import annotations

import math

import numpy as np
from scipy import ndimage

from .detection import robust_location_scale
from .models import SyntheticTrackResult


def shift_and_stack_cutout(
    aligned_frames: list[np.ndarray],
    times_minutes: list[float],
    *,
    center_x: float,
    center_y: float,
    velocity_x_px_min: float,
    velocity_y_px_min: float,
    reference_time_min: float | None = None,
    radius_px: int = 28,
) -> SyntheticTrackResult:
    if len(aligned_frames) != len(times_minutes):
        raise ValueError("Synthetic tracking requires one timestamp per frame.")
    if len(aligned_frames) < 3:
        raise ValueError("Synthetic tracking requires at least three aligned frames.")
    reference_time = (
        float(reference_time_min)
        if reference_time_min is not None
        else float(times_minutes[len(times_minutes) // 2])
    )
    axis = np.arange(-radius_px, radius_px + 1, dtype=np.float64)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    cutouts = []
    skipped = 0
    for frame, time_min in zip(aligned_frames, times_minutes):
        dt = float(time_min) - reference_time
        predicted_x = float(center_x) + float(velocity_x_px_min) * dt
        predicted_y = float(center_y) + float(velocity_y_px_min) * dt
        if (
            predicted_x - radius_px < 0
            or predicted_y - radius_px < 0
            or predicted_x + radius_px >= frame.shape[1]
            or predicted_y + radius_px >= frame.shape[0]
        ):
            skipped += 1
            continue
        sampled = ndimage.map_coordinates(
            np.asarray(frame, dtype=np.float32),
            [predicted_y + yy, predicted_x + xx],
            order=1,
            mode="nearest",
            prefilter=False,
        )
        median, _scale = robust_location_scale(sampled)
        cutouts.append((sampled - median).astype(np.float32))
    if len(cutouts) < 3:
        raise ValueError("Synthetic tracking produced fewer than three usable cutouts.")
    stack = combine_centered_cutouts(cutouts)
    snr, peak_offset = _stack_peak_metrics(stack)
    return SyntheticTrackResult(
        image=stack,
        used_frames=len(cutouts),
        skipped_frames=skipped,
        velocity_x_px_min=float(velocity_x_px_min),
        velocity_y_px_min=float(velocity_y_px_min),
        snr=snr,
        peak_offset_px=peak_offset,
        center_x=float(center_x),
        center_y=float(center_y),
    )


def combine_centered_cutouts(cutouts: list[np.ndarray]) -> np.ndarray:
    if len(cutouts) < 3:
        raise ValueError("At least three centered cutouts are required.")
    cube = np.stack(cutouts, axis=0)
    median = np.median(cube, axis=0)
    deviation = np.abs(cube - median)
    mad = np.median(deviation, axis=0)
    limit = np.maximum(4.0 * 1.4826 * mad, 1e-6)
    valid = deviation <= limit
    stack = np.divide(
        np.sum(np.where(valid, cube, 0.0), axis=0),
        np.maximum(np.sum(valid, axis=0), 1),
    ).astype(np.float32)
    return stack


def velocity_sweep(
    aligned_frames: list[np.ndarray],
    times_minutes: list[float],
    *,
    center_x: float,
    center_y: float,
    velocity_x_values: list[float],
    velocity_y_values: list[float],
    radius_px: int = 20,
    stop_callback=None,
    progress_callback=None,
) -> list[SyntheticTrackResult]:
    results = []
    total = max(1, len(velocity_x_values) * len(velocity_y_values))
    completed = 0
    for vx in velocity_x_values:
        for vy in velocity_y_values:
            if stop_callback and stop_callback():
                return sorted(results, key=lambda item: item.snr, reverse=True)
            result = shift_and_stack_cutout(
                aligned_frames,
                times_minutes,
                center_x=center_x,
                center_y=center_y,
                velocity_x_px_min=vx,
                velocity_y_px_min=vy,
                radius_px=radius_px,
            )
            results.append(result)
            completed += 1
            if progress_callback:
                progress_callback(completed, total)
    return sorted(results, key=lambda item: item.snr, reverse=True)


def _stack_peak_metrics(stack: np.ndarray) -> tuple[float, float]:
    height, width = stack.shape
    center_y, center_x = (height - 1) / 2.0, (width - 1) / 2.0
    yy, xx = np.indices(stack.shape, dtype=np.float32)
    radius = np.hypot(xx - center_x, yy - center_y)
    search = radius <= 4.0
    background = radius >= min(height, width) * 0.28
    median, rms = robust_location_scale(stack[background])
    masked = np.where(search, stack, -np.inf)
    peak_y, peak_x = np.unravel_index(int(np.argmax(masked)), stack.shape)
    peak = float(stack[peak_y, peak_x])
    snr = (peak - median) / max(rms, 1e-6)
    offset = math.hypot(float(peak_x) - center_x, float(peak_y) - center_y)
    return float(snr), float(offset)
