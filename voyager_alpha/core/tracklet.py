from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from .models import Detection, FrameRecord, Tracklet
from .discovery_method import DOCUMENTED_DISCOVERY_METHOD


def link_tracklets(
    detections_by_frame: dict[int, list[Detection]],
    frames: list[FrameRecord],
    max_step_px: float = 35.0,
    min_frames: int = DOCUMENTED_DISCOVERY_METHOD.min_linked_frames,
    min_motion_px_per_frame: float = DOCUMENTED_DISCOVERY_METHOD.min_seed_displacement_px,
    min_median_snr: float = 4.5,
    max_fit_rms_px: float = DOCUMENTED_DISCOVERY_METHOD.borderline_review_rms_px,
    pixel_scale_arcsec: float | None = None,
    match_tolerance_px: float = DOCUMENTED_DISCOVERY_METHOD.match_tolerance_px,
    max_seed_gap_frames: int = 3,
    min_track_occupancy: float = 0.5,
    max_missing_gap_frames: int = 3,
    max_artifact_fraction: float = 0.5,
    strong_fit_rms_px: float = DOCUMENTED_DISCOVERY_METHOD.potential_discovery_rms_px,
    min_motion_arcsec_min: float = 0.017,
    max_motion_arcsec_min: float = 12.5,
) -> list[Tracklet]:
    """Link residual detections with time-aware constant-velocity hypotheses."""

    if not frames or not detections_by_frame:
        return []
    times = _frame_times_minutes(frames)
    valid_indices = sorted(index for index in detections_by_frame if 0 <= index < len(frames))
    cadence = _median_cadence(times)
    if pixel_scale_arcsec is not None and pixel_scale_arcsec > 0:
        min_speed = max(0.0, float(min_motion_arcsec_min)) / pixel_scale_arcsec
        max_speed = max(float(min_motion_arcsec_min), float(max_motion_arcsec_min)) / pixel_scale_arcsec
    else:
        min_speed = float(min_motion_px_per_frame) / max(cadence, 1e-6)
        max_speed = float(max_step_px) / max(cadence, 1e-6)
    tracks_by_signature: dict[tuple[tuple[int, int, int], ...], list[Detection]] = {}

    for position, first_index in enumerate(valid_indices):
        first_detections = detections_by_frame.get(first_index, [])
        if not first_detections:
            continue
        for second_index in valid_indices[position + 1 :]:
            if second_index - first_index > max_seed_gap_frames:
                break
            dt = times[second_index] - times[first_index]
            if dt <= 0:
                continue
            for first in first_detections:
                for second in detections_by_frame.get(second_index, []):
                    vx = (second.x - first.x) / dt
                    vy = (second.y - first.y) / dt
                    speed = math.hypot(vx, vy)
                    if speed < min_speed or speed > max_speed:
                        continue
                    track = _collect_hypothesis(
                        first,
                        vx,
                        vy,
                        times,
                        valid_indices,
                        detections_by_frame,
                        tolerance=match_tolerance_px,
                    )
                    if len(track) < min_frames:
                        continue
                    track, fit = _refine_track(track, times, detections_by_frame, match_tolerance_px)
                    if (
                        len(track) < min_frames
                        or fit[4] > max_fit_rms_px
                        or not _track_has_sequence_support(
                            track,
                            min_occupancy=min_track_occupancy,
                            max_missing_gap_frames=max_missing_gap_frames,
                        )
                        or _maximum_fit_residual(track, times, fit) > max(match_tolerance_px, 1.6 * max_fit_rms_px)
                    ):
                        continue
                    signature = tuple((d.frame_index, round(d.x), round(d.y)) for d in track)
                    current = tracks_by_signature.get(signature)
                    if current is None or _track_quality(track, fit[4]) > _track_quality(current, _fit_motion(current, times)[4]):
                        tracks_by_signature[signature] = track

    candidate_tracks = _remove_overlapping_tracks(list(tracks_by_signature.values()), times)
    tracklets = []
    for track in candidate_tracks:
        x0, y0, vx, vy, fit_rms = _fit_motion(track, times)
        median_snr = float(np.median([d.snr for d in track]))
        if median_snr < min_median_snr or fit_rms > max_fit_rms_px:
            continue
        speed = math.hypot(vx, vy)
        if speed < min_speed:
            continue
        first_time = times[track[0].frame_index]
        last_time = times[track[-1].frame_index]
        arc_minutes = max(0.0, last_time - first_time)
        motion_arcsec = speed * pixel_scale_arcsec if pixel_scale_arcsec else None
        motion_px_per_frame = speed * cadence
        artifact_flags = sorted(
            {
                flag
                for detection in track
                for flag in detection.flags
                if flag in {"star_subtraction_dipole", "near_saturation", "hot_pixel_or_cosmic_ray"}
            }
        )
        artifact_detections = sum(
            bool(set(detection.flags) & {"star_subtraction_dipole", "near_saturation", "hot_pixel_or_cosmic_ray"})
            for detection in track
        )
        if artifact_detections / max(1, len(track)) > max_artifact_fraction:
            continue
        confidence = _confidence(
            len(track),
            fit_rms,
            median_snr,
            len(frames),
            artifact_flags,
            max_fit_rms_px=max_fit_rms_px,
        )
        classification = "unknown_candidate" if fit_rms <= strong_fit_rms_px else "review_candidate"
        position_angle, angle_source = _position_angle(track, times, vx, vy)
        reduced_chi2 = _reduced_chi2(track, times, (x0, y0, vx, vy, fit_rms))
        tracklets.append(
            Tracklet(
                tracklet_id="",
                detections=sorted(track, key=lambda item: item.frame_index),
                frames_detected=len(track),
                motion_px_per_frame=float(motion_px_per_frame),
                motion_arcsec_per_min=float(motion_arcsec) if motion_arcsec is not None else None,
                position_angle_deg=position_angle,
                fit_rms_px=float(fit_rms),
                median_snr=median_snr,
                confidence=confidence,
                velocity_x_px_min=float(vx),
                velocity_y_px_min=float(vy),
                arc_minutes=float(arc_minutes),
                classification=classification,
                artifact_flags=artifact_flags,
                fit_rms_arcsec=float(fit_rms * pixel_scale_arcsec) if pixel_scale_arcsec else None,
                reduced_chi2=reduced_chi2,
                position_angle_source=angle_source,
                speed_regime=_speed_regime(motion_arcsec),
            )
        )
    tracklets.sort(key=lambda item: (item.confidence, item.frames_detected, item.median_snr), reverse=True)
    for index, tracklet in enumerate(tracklets, start=1):
        tracklet.tracklet_id = f"VA-{index:05d}"
    return tracklets


def _collect_hypothesis(
    anchor: Detection,
    vx: float,
    vy: float,
    times: np.ndarray,
    frame_indices: list[int],
    detections_by_frame: dict[int, list[Detection]],
    *,
    tolerance: float,
) -> list[Detection]:
    anchor_time = times[anchor.frame_index]
    track = []
    for frame_index in frame_indices:
        dt = times[frame_index] - anchor_time
        predicted_x = anchor.x + vx * dt
        predicted_y = anchor.y + vy * dt
        candidate = _nearest_detection(
            detections_by_frame.get(frame_index, []),
            predicted_x,
            predicted_y,
            tolerance,
        )
        if candidate is not None:
            track.append(candidate)
    return sorted({d.frame_index: d for d in track}.values(), key=lambda item: item.frame_index)


def _refine_track(
    track: list[Detection],
    times: np.ndarray,
    detections_by_frame: dict[int, list[Detection]],
    tolerance: float,
) -> tuple[list[Detection], tuple[float, float, float, float, float]]:
    current = list(track)
    fit = _fit_motion(current, times)
    for _ in range(3):
        x0, y0, vx, vy, _rms = fit
        origin_time = times[current[0].frame_index]
        rematched = []
        for frame_index in sorted(detections_by_frame):
            dt = times[frame_index] - origin_time
            candidate = _nearest_detection(
                detections_by_frame[frame_index],
                x0 + vx * dt,
                y0 + vy * dt,
                tolerance,
            )
            if candidate is not None:
                rematched.append(candidate)
        rematched = sorted({d.frame_index: d for d in rematched}.values(), key=lambda item: item.frame_index)
        if len(rematched) < 2:
            break
        next_fit = _fit_motion(rematched, times)
        if tuple((d.frame_index, round(d.x, 2), round(d.y, 2)) for d in rematched) == tuple(
            (d.frame_index, round(d.x, 2), round(d.y, 2)) for d in current
        ):
            current, fit = rematched, next_fit
            break
        current, fit = rematched, next_fit
    return current, fit


def _fit_motion(track: list[Detection], times: np.ndarray) -> tuple[float, float, float, float, float]:
    ordered = sorted(track, key=lambda item: item.frame_index)
    sample_times = np.asarray([times[d.frame_index] for d in ordered], dtype=np.float64)
    origin = sample_times[0]
    relative = sample_times - origin
    design = np.column_stack((np.ones(len(relative)), relative))
    xs = np.asarray([d.x for d in ordered], dtype=np.float64)
    ys = np.asarray([d.y for d in ordered], dtype=np.float64)
    sigmas = np.asarray([_centroid_sigma_px(detection) for detection in ordered], dtype=np.float64)
    weights = 1.0 / np.square(sigmas)
    weighted_design = design * np.sqrt(weights)[:, None]
    x_fit = np.linalg.lstsq(weighted_design, xs * np.sqrt(weights), rcond=None)[0]
    y_fit = np.linalg.lstsq(weighted_design, ys * np.sqrt(weights), rcond=None)[0]
    predicted_x = design @ x_fit
    predicted_y = design @ y_fit
    residuals = np.hypot(xs - predicted_x, ys - predicted_y)
    rms = float(np.sqrt(np.mean(np.square(residuals))))
    return float(x_fit[0]), float(y_fit[0]), float(x_fit[1]), float(y_fit[1]), rms


def _maximum_fit_residual(
    track: list[Detection],
    times: np.ndarray,
    fit: tuple[float, float, float, float, float],
) -> float:
    ordered = sorted(track, key=lambda item: item.frame_index)
    x0, y0, vx, vy, _rms = fit
    origin = times[ordered[0].frame_index]
    residuals = [
        math.hypot(
            detection.x - (x0 + vx * (times[detection.frame_index] - origin)),
            detection.y - (y0 + vy * (times[detection.frame_index] - origin)),
        )
        for detection in ordered
    ]
    return max(residuals, default=0.0)


def _track_has_sequence_support(
    track: list[Detection],
    *,
    min_occupancy: float,
    max_missing_gap_frames: int,
) -> bool:
    indices = sorted({detection.frame_index for detection in track})
    if len(indices) < 2:
        return False
    span = indices[-1] - indices[0] + 1
    occupancy = len(indices) / max(1, span)
    if occupancy < float(np.clip(min_occupancy, 0.0, 1.0)):
        return False
    maximum_gap = max((right - left - 1 for left, right in zip(indices, indices[1:])), default=0)
    return maximum_gap <= max(0, int(max_missing_gap_frames))


def _nearest_detection(
    detections: list[Detection],
    x: float,
    y: float,
    tolerance: float,
) -> Detection | None:
    best = None
    best_distance = float(tolerance)
    for detection in detections:
        distance = math.hypot(detection.x - x, detection.y - y)
        if distance <= best_distance:
            best = detection
            best_distance = distance
    return best


def _remove_overlapping_tracks(tracks: list[list[Detection]], times: np.ndarray) -> list[list[Detection]]:
    ranked = sorted(tracks, key=lambda track: _track_quality(track, _fit_motion(track, times)[4]), reverse=True)
    accepted: list[list[Detection]] = []
    used_signatures: list[set[tuple[int, int, int]]] = []
    for track in ranked:
        signature = {(d.frame_index, round(d.x), round(d.y)) for d in track}
        duplicate = False
        for existing in used_signatures:
            overlap = len(signature & existing) / max(1, min(len(signature), len(existing)))
            if overlap >= 0.6:
                duplicate = True
                break
        if not duplicate:
            accepted.append(track)
            used_signatures.append(signature)
    return accepted


def _track_quality(track: list[Detection], fit_rms: float) -> float:
    return len(track) * 10.0 + float(np.median([d.snr for d in track])) - fit_rms * 4.0


def _confidence(
    frames_detected: int,
    fit_rms_px: float,
    median_snr: float,
    total_frames: int,
    artifact_flags: list[str],
    *,
    max_fit_rms_px: float,
) -> float:
    required_span = min(5, max(3, total_frames))
    frame_score = min(1.0, frames_detected / required_span)
    fit_score = max(0.0, 1.0 - fit_rms_px / max(float(max_fit_rms_px), 1e-6))
    snr_score = min(1.0, median_snr / 12.0)
    penalty = min(0.25, 0.08 * len(artifact_flags))
    return round(max(0.0, min(0.99, 0.45 * frame_score + 0.35 * fit_score + 0.20 * snr_score - penalty)), 3)


def _frame_times_minutes(frames: list[FrameRecord]) -> np.ndarray:
    if any(frame.midpoint_jd is None or not np.isfinite(frame.midpoint_jd) for frame in frames):
        raise ValueError("Tracklet linking requires a valid exposure midpoint for every frame.")
    base = float(frames[0].midpoint_jd)
    values = [(float(frame.midpoint_jd) - base) * 1440.0 for frame in frames]
    result = np.asarray(values, dtype=np.float64)
    if len(result) > 1 and np.any(np.diff(result) <= 0):
        raise ValueError("Tracklet linking requires strictly increasing exposure times.")
    return result


def _centroid_sigma_px(detection: Detection) -> float:
    if detection.snr <= 0 or detection.fwhm_px <= 0:
        return 1.0
    return max(0.05, float(detection.fwhm_px) / (2.355 * float(detection.snr)))


def _reduced_chi2(
    track: list[Detection],
    times: np.ndarray,
    fit: tuple[float, float, float, float, float],
) -> float | None:
    if len(track) <= 2:
        return None
    ordered = sorted(track, key=lambda item: item.frame_index)
    x0, y0, vx, vy, _rms = fit
    origin = times[ordered[0].frame_index]
    chi2 = 0.0
    for detection in ordered:
        dt = times[detection.frame_index] - origin
        residual = math.hypot(detection.x - (x0 + vx * dt), detection.y - (y0 + vy * dt))
        chi2 += (residual / _centroid_sigma_px(detection)) ** 2
    return float(chi2 / max(1, 2 * len(ordered) - 4))


def _position_angle(
    track: list[Detection],
    times: np.ndarray,
    fallback_vx: float,
    fallback_vy: float,
) -> tuple[float, str]:
    ordered = sorted(track, key=lambda item: item.frame_index)
    if len(ordered) >= 2 and all(item.ra is not None and item.dec is not None for item in ordered):
        sample_times = np.asarray([times[item.frame_index] for item in ordered], dtype=np.float64)
        relative = sample_times - sample_times[0]
        design = np.column_stack((np.ones(len(relative)), relative))
        ra = np.rad2deg(np.unwrap(np.deg2rad([float(item.ra) for item in ordered])))
        dec = np.asarray([float(item.dec) for item in ordered], dtype=np.float64)
        east = (ra - ra[0]) * math.cos(math.radians(float(np.mean(dec)))) * 3600.0
        north = (dec - dec[0]) * 3600.0
        east_rate = float(np.linalg.lstsq(design, east, rcond=None)[0][1])
        north_rate = float(np.linalg.lstsq(design, north, rcond=None)[0][1])
        if math.hypot(east_rate, north_rate) > 0:
            return float((math.degrees(math.atan2(east_rate, north_rate)) + 360.0) % 360.0), "sky_wcs"
    return float((math.degrees(math.atan2(fallback_vx, fallback_vy)) + 360.0) % 360.0), "pixel_grid"


def _speed_regime(motion_arcsec_min: float | None) -> str:
    if motion_arcsec_min is None:
        return "pixel_only"
    if motion_arcsec_min < 0.1:
        return "very_slow"
    if motion_arcsec_min <= 5.0:
        return "nominal"
    return "fast"


def _median_cadence(times: np.ndarray) -> float:
    if len(times) < 2:
        return 1.0
    differences = np.diff(times)
    differences = differences[differences > 0]
    return float(np.median(differences)) if differences.size else 1.0


def tracklets_by_frame(tracklets: list[Tracklet]) -> dict[int, list[Tracklet]]:
    result = defaultdict(list)
    for tracklet in tracklets:
        for detection in tracklet.detections:
            result[detection.frame_index].append(tracklet)
    return dict(result)
