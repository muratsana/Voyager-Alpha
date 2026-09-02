from __future__ import annotations

import numpy as np
import sep
from scipy import ndimage

from .models import RegistrationSolution


def normalize_for_registration(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    lo, hi = np.nanpercentile(data, (1.0, 99.8))
    if not np.isfinite(hi) or hi <= lo:
        lo = float(np.nanmedian(data))
        hi = float(np.nanmax(data))
    if not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(data, dtype=np.float32)
    clipped = np.clip(data, lo, hi)
    clipped -= np.nanmedian(clipped)
    scale = 1.4826 * np.nanmedian(np.abs(clipped - np.nanmedian(clipped)))
    if not np.isfinite(scale) or scale <= 0:
        scale = np.nanstd(clipped)
    if not np.isfinite(scale) or scale <= 0:
        return np.zeros_like(clipped, dtype=np.float32)
    return np.nan_to_num(clipped / scale, copy=False).astype(np.float32)


def _registration_view(image: np.ndarray, max_dim: int = 1400) -> tuple[np.ndarray, int]:
    data = np.asarray(image, dtype=np.float32)
    factor = max(1, int(np.ceil(max(data.shape) / max_dim)))
    return data[::factor, ::factor], factor


def estimate_subpixel_shift(
    reference: np.ndarray,
    moving: np.ndarray,
    max_shift_px: float = 180.0,
) -> tuple[float, float, float]:
    """Estimate the translation that moves ``moving`` onto ``reference``."""

    if reference.shape != moving.shape:
        raise ValueError("Registration requires frames with identical shapes.")
    ref_view, factor = _registration_view(reference)
    mov_view, _ = _registration_view(moving)
    ref = normalize_for_registration(ref_view)
    mov = normalize_for_registration(mov_view)
    cross_power = np.fft.fft2(ref) * np.fft.fft2(mov).conj()
    denom = np.abs(cross_power)
    cross_power = np.divide(cross_power, denom, out=np.zeros_like(cross_power), where=denom > 0)
    corr = np.abs(np.fft.ifft2(cross_power))
    peak_y, peak_x = np.unravel_index(int(np.argmax(corr)), corr.shape)
    height, width = corr.shape
    if peak_y > height // 2:
        peak_y -= height
    if peak_x > width // 2:
        peak_x -= width
    peak_y = float(peak_y) + _parabolic_offset(corr[:, int(peak_x) % width], int(peak_y) % height)
    peak_x = float(peak_x) + _parabolic_offset(corr[int(peak_y) % height, :], int(peak_x) % width)
    dy = float(np.clip(peak_y * factor, -max_shift_px, max_shift_px))
    dx = float(np.clip(peak_x * factor, -max_shift_px, max_shift_px))
    return dy, dx, float(np.max(corr))


def _parabolic_offset(values: np.ndarray, index: int) -> float:
    length = len(values)
    left = float(values[(index - 1) % length])
    center = float(values[index % length])
    right = float(values[(index + 1) % length])
    denominator = left - 2.0 * center + right
    if abs(denominator) < 1e-12:
        return 0.0
    return float(np.clip(0.5 * (left - right) / denominator, -0.75, 0.75))


def estimate_integer_shift(
    reference: np.ndarray,
    moving: np.ndarray,
    max_shift_px: int = 180,
) -> tuple[int, int, float]:
    dy, dx, peak = estimate_subpixel_shift(reference, moving, max_shift_px=float(max_shift_px))
    return int(round(dy)), int(round(dx)), peak


def detect_registration_stars(
    image: np.ndarray,
    *,
    threshold_sigma: float = 7.0,
    max_stars: int = 350,
    edge_margin: int = 12,
) -> np.ndarray:
    data = _native_float32(image)
    background = sep.Background(data)
    residual = data - background
    objects = sep.extract(
        residual,
        threshold_sigma,
        err=max(float(background.globalrms), 1e-6),
        minarea=5,
    )
    if len(objects) == 0:
        return np.empty((0, 2), dtype=np.float64)
    height, width = data.shape
    valid = (
        (objects["x"] >= edge_margin)
        & (objects["x"] < width - edge_margin)
        & (objects["y"] >= edge_margin)
        & (objects["y"] < height - edge_margin)
        & (objects["a"] > 0.5)
        & (objects["b"] > 0.5)
        & (objects["a"] / np.maximum(objects["b"], 1e-6) < 2.4)
    )
    objects = objects[valid]
    if len(objects) == 0:
        return np.empty((0, 2), dtype=np.float64)
    order = np.argsort(objects["flux"])[::-1][:max_stars]
    return np.column_stack((objects["x"][order], objects["y"][order])).astype(np.float64)


def register_frame(
    reference: np.ndarray,
    moving: np.ndarray,
    *,
    max_shift_px: float = 180.0,
    match_radius_px: float = 7.0,
) -> tuple[np.ndarray, RegistrationSolution]:
    if reference.shape != moving.shape:
        raise ValueError("Registration requires frames with identical shapes.")

    ref_stars = detect_registration_stars(reference)
    mov_stars = detect_registration_stars(moving)
    candidates: list[RegistrationSolution] = []
    for orientation, oriented, oriented_stars, orientation_matrix, orientation_offset in _orientation_hypotheses(
        moving,
        mov_stars,
    ):
        refinement = _solve_oriented_registration(
            reference,
            oriented,
            ref_stars,
            oriented_stars,
            max_shift_px=max_shift_px,
            match_radius_px=match_radius_px,
        )
        if orientation != "identity" and not refinement.method.startswith("star"):
            continue
        solution = _compose_registration(
            refinement,
            orientation_matrix,
            orientation_offset,
            orientation,
        )
        if refinement.method.startswith("star") and not _is_plausible_refinement(refinement.matrix_xy):
            continue
        candidates.append(solution)

    if not candidates:
        raise ValueError("No valid registration orientation could be evaluated.")
    solution = max(candidates, key=_registration_candidate_score)
    fill_value = float(np.nanmedian(reference))
    aligned = warp_affine(moving, solution, reference.shape, fill_value=fill_value)
    return _refine_star_residual(reference, moving, aligned, solution, ref_stars, fill_value)


def _orientation_hypotheses(
    moving: np.ndarray,
    moving_stars: np.ndarray,
) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    height, width = moving.shape
    rotated_stars = np.asarray([width - 1.0, height - 1.0], dtype=np.float64) - moving_stars
    return [
        (
            "identity",
            np.asarray(moving, dtype=np.float32),
            moving_stars,
            np.eye(2, dtype=np.float64),
            np.zeros(2, dtype=np.float64),
        ),
        (
            "rot180",
            np.ascontiguousarray(np.rot90(moving, 2), dtype=np.float32),
            rotated_stars,
            -np.eye(2, dtype=np.float64),
            np.asarray([width - 1.0, height - 1.0], dtype=np.float64),
        ),
    ]


def _solve_oriented_registration(
    reference: np.ndarray,
    moving: np.ndarray,
    ref_stars: np.ndarray,
    mov_stars: np.ndarray,
    *,
    max_shift_px: float,
    match_radius_px: float,
) -> RegistrationSolution:
    dy, dx, peak = estimate_subpixel_shift(reference, moving, max_shift_px=max_shift_px)
    pairs = _match_star_pairs(ref_stars, mov_stars, dx, dy, match_radius_px)
    if len(pairs) >= 8:
        moving_points = np.asarray([item[0] for item in pairs], dtype=np.float64)
        reference_points = np.asarray([item[1] for item in pairs], dtype=np.float64)
        matrix, offset, rms, kept = _robust_affine_fit(moving_points, reference_points)
        solution = RegistrationSolution(
            matrix_xy=matrix,
            offset_xy=offset,
            rms_px=rms,
            matched_stars=kept,
            method="star-affine",
            phase_peak=peak,
        )
    elif len(pairs) >= 3:
        offsets = np.asarray([reference_point - moving_point for moving_point, reference_point in pairs])
        offset, residuals, keep = _robust_translation(offsets)
        solution = RegistrationSolution(
            matrix_xy=np.eye(2, dtype=np.float64),
            offset_xy=offset,
            rms_px=float(np.sqrt(np.mean(np.square(residuals)))),
            matched_stars=int(np.count_nonzero(keep)),
            method="star-translation",
            phase_peak=peak,
        )
    else:
        solution = RegistrationSolution(
            matrix_xy=np.eye(2, dtype=np.float64),
            offset_xy=np.asarray([dx, dy], dtype=np.float64),
            rms_px=2.5,
            matched_stars=len(pairs),
            method="phase-subpixel",
            phase_peak=peak,
        )
    return solution


def _robust_translation(offsets: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(offsets, dtype=np.float64)
    offset = np.median(values, axis=0)
    residuals = np.hypot(*((values - offset).T))
    residual_median = float(np.median(residuals))
    residual_mad = float(np.median(np.abs(residuals - residual_median)))
    limit = max(1.25, residual_median + 3.5 * max(1.4826 * residual_mad, 0.1))
    keep = residuals <= limit
    if np.count_nonzero(keep) >= 3:
        offset = np.median(values[keep], axis=0)
        residuals = np.hypot(*((values[keep] - offset).T))
    else:
        keep = np.ones(len(values), dtype=bool)
    return offset, residuals, keep


def _compose_registration(
    refinement: RegistrationSolution,
    orientation_matrix: np.ndarray,
    orientation_offset: np.ndarray,
    orientation: str,
) -> RegistrationSolution:
    matrix = refinement.matrix_xy @ orientation_matrix
    offset = refinement.matrix_xy @ orientation_offset + refinement.offset_xy
    method = refinement.method if orientation == "identity" else f"{refinement.method}-{orientation}"
    return RegistrationSolution(
        matrix_xy=matrix,
        offset_xy=offset,
        rms_px=refinement.rms_px,
        matched_stars=refinement.matched_stars,
        method=method,
        phase_peak=refinement.phase_peak,
    )


def _is_plausible_refinement(matrix: np.ndarray) -> bool:
    singular_values = np.linalg.svd(np.asarray(matrix, dtype=np.float64), compute_uv=False)
    determinant = float(np.linalg.det(matrix))
    return bool(
        np.all(np.isfinite(singular_values))
        and 0.85 <= float(np.min(singular_values))
        and float(np.max(singular_values)) <= 1.15
        and 0.75 <= determinant <= 1.25
    )


def _registration_candidate_score(solution: RegistrationSolution) -> tuple[float, ...]:
    is_star_solution = float(solution.method.startswith("star"))
    identity_tie_break = float("rot180" not in solution.method)
    return (
        is_star_solution,
        float(solution.matched_stars),
        -float(solution.rms_px),
        float(solution.phase_peak),
        identity_tie_break,
    )


def _refine_star_residual(
    reference: np.ndarray,
    moving: np.ndarray,
    aligned: np.ndarray,
    solution: RegistrationSolution,
    ref_stars: np.ndarray,
    fill_value: float,
) -> tuple[np.ndarray, RegistrationSolution]:
    if not solution.method.startswith("star") or len(ref_stars) < 3:
        return aligned, solution
    aligned_stars = detect_registration_stars(aligned)
    pairs = _match_star_pairs(ref_stars, aligned_stars, 0.0, 0.0, 2.5)
    if len(pairs) < 3:
        return aligned, solution
    offsets = np.asarray([reference_point - moving_point for moving_point, reference_point in pairs])
    correction, _residuals, keep = _robust_translation(offsets)
    if np.count_nonzero(keep) < 3 or float(np.hypot(*correction)) > 2.5:
        return aligned, solution

    refined = RegistrationSolution(
        matrix_xy=solution.matrix_xy,
        offset_xy=solution.offset_xy + correction,
        rms_px=solution.rms_px,
        matched_stars=solution.matched_stars,
        method=solution.method,
        phase_peak=solution.phase_peak,
    )
    aligned = warp_affine(moving, refined, reference.shape, fill_value=fill_value)
    final_stars = detect_registration_stars(aligned)
    final_pairs = _match_star_pairs(ref_stars, final_stars, 0.0, 0.0, 2.0)
    if len(final_pairs) >= 3:
        final_offsets = np.asarray(
            [reference_point - moving_point for moving_point, reference_point in final_pairs]
        )
        _offset, residuals, final_keep = _robust_translation(final_offsets)
        refined.rms_px = float(np.sqrt(np.mean(np.square(residuals))))
        refined.matched_stars = max(refined.matched_stars, int(np.count_nonzero(final_keep)))
    return aligned, refined


def _match_star_pairs(
    reference_stars: np.ndarray,
    moving_stars: np.ndarray,
    dx: float,
    dy: float,
    radius: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    if len(reference_stars) == 0 or len(moving_stars) == 0:
        return []
    shifted = moving_stars + np.asarray([dx, dy], dtype=np.float64)
    delta = shifted[:, None, :] - reference_stars[None, :, :]
    distances = np.hypot(delta[:, :, 0], delta[:, :, 1])
    candidates = []
    for moving_index in range(len(moving_stars)):
        reference_index = int(np.argmin(distances[moving_index]))
        distance = float(distances[moving_index, reference_index])
        if distance <= radius:
            candidates.append((distance, moving_index, reference_index))
    pairs = []
    used_moving = set()
    used_reference = set()
    for _distance, moving_index, reference_index in sorted(candidates):
        if moving_index in used_moving or reference_index in used_reference:
            continue
        used_moving.add(moving_index)
        used_reference.add(reference_index)
        pairs.append((moving_stars[moving_index], reference_stars[reference_index]))
    return pairs


def _robust_affine_fit(
    moving_points: np.ndarray,
    reference_points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, int]:
    keep = np.ones(len(moving_points), dtype=bool)
    matrix = np.eye(2, dtype=np.float64)
    offset = np.zeros(2, dtype=np.float64)
    for _ in range(5):
        source = moving_points[keep]
        target = reference_points[keep]
        if len(source) < 6:
            break
        design = np.column_stack((source, np.ones(len(source))))
        coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
        matrix = coefficients[:2, :].T
        offset = coefficients[2, :]
        predicted = moving_points @ matrix.T + offset
        residuals = np.hypot(*(predicted - reference_points).T)
        median = float(np.median(residuals[keep]))
        mad = float(np.median(np.abs(residuals[keep] - median)))
        limit = max(1.25, median + 3.5 * max(1.4826 * mad, 0.1))
        next_keep = residuals <= min(limit, 4.0)
        if np.array_equal(next_keep, keep):
            break
        keep = next_keep
    predicted = moving_points @ matrix.T + offset
    residuals = np.hypot(*(predicted - reference_points).T)
    rms = float(np.sqrt(np.mean(np.square(residuals[keep])))) if np.any(keep) else float("inf")
    return matrix, offset, rms, int(np.count_nonzero(keep))


def warp_affine(
    image: np.ndarray,
    solution: RegistrationSolution,
    output_shape: tuple[int, int],
    *,
    fill_value: float = 0.0,
) -> np.ndarray:
    forward = np.asarray(solution.matrix_xy, dtype=np.float64)
    inverse = np.linalg.inv(forward)
    inverse_offset = -inverse @ np.asarray(solution.offset_xy, dtype=np.float64)
    swap = np.asarray([[0.0, 1.0], [1.0, 0.0]])
    matrix_yx = swap @ inverse @ swap
    offset_yx = swap @ inverse_offset
    return ndimage.affine_transform(
        np.asarray(image, dtype=np.float32),
        matrix_yx,
        offset=offset_yx,
        output_shape=output_shape,
        order=1,
        mode="constant",
        cval=float(fill_value),
        prefilter=False,
    ).astype(np.float32, copy=False)


def shift_image_integer(image: np.ndarray, dy: int, dx: int, fill_value: float = 0.0) -> np.ndarray:
    solution = RegistrationSolution(
        matrix_xy=np.eye(2, dtype=np.float64),
        offset_xy=np.asarray([float(dx), float(dy)]),
        rms_px=0.0,
        matched_stars=0,
        method="integer",
    )
    return warp_affine(image, solution, np.asarray(image).shape, fill_value=fill_value)


def align_frames(
    frames: list[np.ndarray],
    progress_callback=None,
) -> tuple[list[np.ndarray], list[dict[str, float | int | str]]]:
    if not frames:
        return [], []
    reference = np.asarray(frames[0], dtype=np.float32)
    aligned = [reference]
    metrics = [{"dy": 0.0, "dx": 0.0, "peak": 1.0, "rms": 0.0, "matched": 0, "method": "reference"}]
    for index, frame in enumerate(frames[1:], start=2):
        if progress_callback:
            progress_callback(index, len(frames))
        registered, solution = register_frame(reference, frame)
        aligned.append(registered)
        metrics.append(
            {
                "dy": float(solution.offset_xy[1]),
                "dx": float(solution.offset_xy[0]),
                "peak": solution.phase_peak,
                "rms": solution.rms_px,
                "matched": solution.matched_stars,
                "method": solution.method,
            }
        )
    return aligned, metrics


def _native_float32(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    if not data.dtype.isnative:
        data = data.byteswap().view(data.dtype.newbyteorder("="))
    return np.ascontiguousarray(data)
