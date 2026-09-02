from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

try:
    from exoplanet_core import quad_limbdark_light_curve
except ImportError:  # pragma: no cover - diagnostics and packaging require it
    quad_limbdark_light_curve = None


@dataclass
class ApertureMeasurement:
    flux: float
    uncertainty: float
    x: float
    y: float
    background: float
    background_rms: float
    aperture_pixels: int
    flags: list[str] = field(default_factory=list)
    peak_adu: float = float("nan")
    linearity_limit_adu: float | None = None
    saturated_pixels: int = 0
    fwhm_px: float = float("nan")
    snr: float = float("nan")


@dataclass
class TransitModelFit:
    success: bool
    model_flux: list[float]
    mid_transit_jd: float | None = None
    duration_minutes: float | None = None
    radius_ratio: float | None = None
    impact_parameter: float | None = None
    depth: float | None = None
    snr: float | None = None
    reduced_chi2: float | None = None
    delta_bic: float | None = None
    message: str = ""


@dataclass
class LightCurveResult:
    times_jd: list[float]
    relative_flux: list[float]
    target_flux: list[float]
    comparison_flux: list[float]
    scatter: float
    depth: float
    transit_candidate: bool
    message: str
    raw_relative_flux: list[float] = field(default_factory=list)
    flux_uncertainty: list[float] = field(default_factory=list)
    valid_mask: list[bool] = field(default_factory=list)
    comparison_weights: list[float] = field(default_factory=list)
    quality_flags: list[str] = field(default_factory=list)
    depth_uncertainty: float = float("nan")
    detection_snr: float = 0.0
    duration_minutes: float | None = None
    mid_transit_jd: float | None = None
    model_fit: TransitModelFit | None = None
    time_system: str = "JD_UTC"


def recenter_star(
    image: np.ndarray,
    x: float,
    y: float,
    *,
    search_radius: int = 5,
    max_offset: float = 4.0,
) -> tuple[float, float, list[str]]:
    data = np.asarray(image, dtype=np.float32)
    x0 = max(0, int(np.floor(x)) - search_radius)
    x1 = min(data.shape[1], int(np.floor(x)) + search_radius + 1)
    y0 = max(0, int(np.floor(y)) - search_radius)
    y1 = min(data.shape[0], int(np.floor(y)) + search_radius + 1)
    flags = []
    if x1 - x0 < 5 or y1 - y0 < 5:
        return float(x), float(y), ["centroid_edge"]
    patch = data[y0:y1, x0:x1]
    finite = patch[np.isfinite(patch)]
    if finite.size < 12:
        return float(x), float(y), ["centroid_invalid"]
    background = float(np.nanpercentile(finite, 25.0))
    weights = np.clip(patch - background, 0.0, None)
    peak = float(np.nanmax(weights))
    if not np.isfinite(peak) or peak <= 0:
        return float(x), float(y), ["centroid_no_signal"]
    weights[weights < peak * 0.08] = 0.0
    total = float(np.nansum(weights))
    if total <= 0:
        return float(x), float(y), ["centroid_no_signal"]
    yy, xx = np.indices(patch.shape, dtype=np.float64)
    measured_x = float(np.nansum((xx + x0) * weights) / total)
    measured_y = float(np.nansum((yy + y0) * weights) / total)
    offset = float(np.hypot(measured_x - x, measured_y - y))
    if offset > max_offset:
        flags.append("centroid_large_shift")
        return float(x), float(y), flags
    return measured_x, measured_y, flags


def aperture_measurement(
    image: np.ndarray,
    x: float,
    y: float,
    *,
    aperture_radius: float = 6.0,
    annulus_inner: float | None = None,
    annulus_outer: float | None = None,
    recenter: bool = True,
    linearity_limit_adu: float | None = None,
    gain_e_per_adu: float | None = None,
) -> ApertureMeasurement:
    data = np.asarray(image, dtype=np.float32)
    annulus_inner = float(annulus_inner or max(aperture_radius + 3.0, aperture_radius * 1.5))
    annulus_outer = float(annulus_outer or max(annulus_inner + 4.0, aperture_radius * 2.3))
    flags = []
    measured_x, measured_y = float(x), float(y)
    if recenter:
        measured_x, measured_y, centroid_flags = recenter_star(data, x, y)
        flags.extend(centroid_flags)
    y0 = max(0, int(np.floor(measured_y - annulus_outer - 1)))
    y1 = min(data.shape[0], int(np.ceil(measured_y + annulus_outer + 2)))
    x0 = max(0, int(np.floor(measured_x - annulus_outer - 1)))
    x1 = min(data.shape[1], int(np.ceil(measured_x + annulus_outer + 2)))
    if x1 <= x0 or y1 <= y0:
        return ApertureMeasurement(float("nan"), float("nan"), measured_x, measured_y, 0.0, 0.0, 0, [*flags, "outside_frame"])
    patch = data[y0:y1, x0:x1]
    yy, xx = np.indices(patch.shape, dtype=np.float32)
    radius = np.hypot(xx + x0 - measured_x, yy + y0 - measured_y)
    aperture = radius <= aperture_radius
    annulus = (radius >= annulus_inner) & (radius <= annulus_outer)
    aperture_pixels = int(np.count_nonzero(aperture))
    annulus_values = patch[annulus & np.isfinite(patch)]
    if aperture_pixels < 4 or annulus_values.size < 12:
        return ApertureMeasurement(float("nan"), float("nan"), measured_x, measured_y, 0.0, 0.0, aperture_pixels, [*flags, "insufficient_pixels"])
    background = float(np.nanmedian(annulus_values))
    mad = float(np.nanmedian(np.abs(annulus_values - background)))
    background_rms = max(1.4826 * mad, 1e-6)
    source_values = patch[aperture]
    flux = float(np.nansum(source_values - background))
    gain = float(gain_e_per_adu) if gain_e_per_adu and gain_e_per_adu > 0 else 1.0
    background_variance = aperture_pixels * background_rms**2 * (
        1.0 + aperture_pixels / max(int(annulus_values.size), 1)
    )
    uncertainty = float(np.sqrt(max(flux, 0.0) / gain + background_variance))
    if not np.isfinite(flux) or flux <= 0:
        flags.append("nonpositive_flux")
    peak_adu = float(np.nanmax(source_values))
    saturated_pixels = 0
    if linearity_limit_adu is not None and np.isfinite(linearity_limit_adu) and linearity_limit_adu > 0:
        saturated_pixels = int(np.count_nonzero(source_values >= float(linearity_limit_adu)))
        if saturated_pixels:
            flags.append("saturated")
        elif peak_adu >= float(linearity_limit_adu) * 0.90:
            flags.append("near_linearity_limit")
    weights = np.clip(patch - background, 0.0, None) * aperture
    weight_sum = float(np.nansum(weights))
    if weight_sum > 0:
        variance_x = float(np.nansum(np.square(xx + x0 - measured_x) * weights) / weight_sum)
        variance_y = float(np.nansum(np.square(yy + y0 - measured_y) * weights) / weight_sum)
        fwhm_px = 2.35482 * np.sqrt(max((variance_x + variance_y) / 2.0, 0.0))
    else:
        fwhm_px = float("nan")
    snr = flux / uncertainty if uncertainty > 0 else float("nan")
    return ApertureMeasurement(
        flux=flux,
        uncertainty=uncertainty,
        x=measured_x,
        y=measured_y,
        background=background,
        background_rms=background_rms,
        aperture_pixels=aperture_pixels,
        flags=flags,
        peak_adu=peak_adu,
        linearity_limit_adu=linearity_limit_adu,
        saturated_pixels=saturated_pixels,
        fwhm_px=float(fwhm_px),
        snr=float(snr),
    )


def aperture_flux(
    image: np.ndarray,
    x: float,
    y: float,
    *,
    aperture_radius: float = 6.0,
    annulus_inner: float = 9.0,
    annulus_outer: float = 14.0,
) -> float:
    return aperture_measurement(
        image,
        x,
        y,
        aperture_radius=aperture_radius,
        annulus_inner=annulus_inner,
        annulus_outer=annulus_outer,
        recenter=False,
    ).flux


def differential_light_curve(
    frames: list[np.ndarray],
    times_jd: list[float],
    target_xy: tuple[float, float],
    comparison_xy: list[tuple[float, float]],
    *,
    aperture_radius: float = 6.0,
    detrend_order: int = 2,
) -> LightCurveResult:
    if len(frames) < 5:
        raise ValueError("Transit fotometrisi için en az 5 kare gerekir.")
    if len(times_jd) != len(frames):
        raise ValueError("Her kare için geçerli zaman bilgisi gerekir.")
    if not comparison_xy:
        raise ValueError("En az bir karşılaştırma yıldızı seçilmelidir.")

    target_values = []
    target_errors = []
    comparison_values = []
    comparison_errors = []
    for frame in frames:
        target = aperture_measurement(frame, *target_xy, aperture_radius=aperture_radius)
        comparisons = [
            aperture_measurement(frame, *xy, aperture_radius=aperture_radius)
            for xy in comparison_xy
        ]
        target_values.append(target.flux)
        target_errors.append(target.uncertainty)
        comparison_values.append([value.flux for value in comparisons])
        comparison_errors.append([value.uncertainty for value in comparisons])

    return differential_light_curve_from_fluxes(
        times_jd,
        target_values,
        comparison_values,
        target_uncertainties=target_errors,
        comparison_uncertainties=comparison_errors,
        detrend_order=detrend_order,
    )


def differential_light_curve_from_fluxes(
    times_jd: list[float],
    target_values: list[float],
    comparison_values,
    *,
    target_uncertainties=None,
    comparison_uncertainties=None,
    detrend_order: int = 2,
    time_system: str = "JD_UTC",
) -> LightCurveResult:
    times = np.asarray(times_jd, dtype=np.float64)
    target = np.asarray(target_values, dtype=np.float64)
    comparisons = np.asarray(comparison_values, dtype=np.float64)
    if comparisons.ndim == 1:
        comparisons = comparisons[:, None]
    if len(times) != len(target) or comparisons.shape[0] != len(times):
        raise ValueError("Zaman, hedef ve karşılaştırma ölçüm sayıları eşit olmalıdır.")
    if len(times) < 5:
        raise ValueError("Transit fotometrisi için en az 5 ölçüm gerekir.")
    if comparisons.shape[1] == 0:
        raise ValueError("En az bir karşılaştırma yıldızı gerekir.")

    target_error = _coerce_uncertainties(target_uncertainties, target)
    comparison_error = _coerce_comparison_uncertainties(comparison_uncertainties, comparisons)
    ensemble, ensemble_error, weights, comparison_flags = _weighted_comparison_ensemble(
        comparisons,
        comparison_error,
    )
    target_median = _positive_median(target)
    target_normalized = target / target_median
    raw_relative = target_normalized / ensemble
    relative_error = np.abs(raw_relative) * np.sqrt(
        np.square(target_error / np.maximum(np.abs(target), 1e-12))
        + np.square(ensemble_error / np.maximum(np.abs(ensemble), 1e-12))
    )
    valid = (
        np.isfinite(times)
        & np.isfinite(raw_relative)
        & np.isfinite(relative_error)
        & (raw_relative > 0)
        & (relative_error > 0)
    )
    if int(np.count_nonzero(valid)) < 5:
        raise ValueError("Geçerli diferansiyel fotometri ölçümü için yeterli kare yok.")

    preliminary_flux = raw_relative / float(np.nanmedian(raw_relative[valid]))
    preliminary_box = _search_single_transit(times, preliminary_flux, relative_error, valid)
    excluded_from_trend = np.zeros(len(times), dtype=bool)
    if preliminary_box["snr"] >= 3.0:
        excluded_from_trend[np.asarray(preliminary_box["in_transit_indices"], dtype=int)] = True
    trend = _robust_polynomial_trend(
        times,
        raw_relative,
        relative_error,
        valid,
        detrend_order,
        exclude=excluded_from_trend,
    )
    relative = raw_relative / trend
    relative_error = relative_error / np.maximum(np.abs(trend), 1e-12)
    baseline = float(np.nanmedian(relative[valid]))
    relative /= baseline
    relative_error /= baseline
    scatter = _robust_sigma(relative[valid] - 1.0)
    box = _search_single_transit(times, relative, relative_error, valid)
    model_fit = fit_limb_darkened_transit(
        times,
        relative,
        relative_error,
        valid,
        initial_midpoint=box["midpoint_jd"],
        initial_duration_minutes=box["duration_minutes"],
        initial_depth=box["depth"],
    )

    depth = float(max(box["depth"], model_fit.depth or 0.0 if model_fit.success else box["depth"]))
    detection_snr = float(max(box["snr"], model_fit.snr or 0.0 if model_fit.success else box["snr"]))
    coverage_ok = bool(box["coverage_ok"])
    model_supported = bool(model_fit.success and (model_fit.delta_bic or 0.0) >= 6.0)
    candidate = bool(
        coverage_ok
        and box["in_transit_points"] >= 2
        and depth >= max(0.002, 3.0 * box["depth_uncertainty"])
        and detection_snr >= 5.0
        and (model_supported or box["snr"] >= 7.0)
    )
    normalized_time_system = str(time_system or "JD_UTC").upper()
    quality_flags = list(comparison_flags)
    if normalized_time_system != "BJD_TDB":
        quality_flags.append("time_is_jd_utc_not_bjd_tdb")
    if not coverage_ok:
        quality_flags.append("incomplete_transit_baseline")
    if scatter > 0.02:
        quality_flags.append("high_photometric_scatter")
    if not model_fit.success:
        quality_flags.append("physical_fit_unavailable")
    message = (
        f"Transit adayı: derinlik {depth * 100:.3f}%, SNR {detection_snr:.1f}, "
        f"süre {box['duration_minutes']:.1f} dk."
        if candidate
        else f"Belirgin transit yok: en güçlü düşüş {depth * 100:.3f}%, SNR {detection_snr:.1f}."
    )
    return LightCurveResult(
        times_jd=[float(value) for value in times],
        relative_flux=[float(value) for value in relative],
        target_flux=[float(value) for value in target],
        comparison_flux=[float(value) for value in ensemble],
        scatter=float(scatter),
        depth=depth,
        transit_candidate=candidate,
        message=message,
        raw_relative_flux=[float(value) for value in raw_relative],
        flux_uncertainty=[float(value) for value in relative_error],
        valid_mask=[bool(value) for value in valid],
        comparison_weights=[float(value) for value in weights],
        quality_flags=quality_flags,
        depth_uncertainty=float(box["depth_uncertainty"]),
        detection_snr=detection_snr,
        duration_minutes=float(box["duration_minutes"]),
        mid_transit_jd=float(box["midpoint_jd"]),
        model_fit=model_fit,
        time_system=normalized_time_system,
    )


def fit_limb_darkened_transit(
    times_jd,
    relative_flux,
    uncertainties,
    valid_mask=None,
    *,
    initial_midpoint=None,
    initial_duration_minutes=None,
    initial_depth=None,
    limb_darkening=(0.3, 0.2),
) -> TransitModelFit:
    times = np.asarray(times_jd, dtype=np.float64)
    flux = np.asarray(relative_flux, dtype=np.float64)
    error = np.asarray(uncertainties, dtype=np.float64)
    valid = np.asarray(valid_mask if valid_mask is not None else np.ones(len(times), dtype=bool), dtype=bool)
    valid &= np.isfinite(times) & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    empty_model = np.full(len(times), np.nan, dtype=np.float64)
    if quad_limbdark_light_curve is None:
        return TransitModelFit(False, empty_model.tolist(), message="exoplanet-core is not available")
    if int(np.count_nonzero(valid)) < 6:
        return TransitModelFit(False, empty_model.tolist(), message="Not enough valid points for physical fit")
    t_ref = float(np.nanmedian(times[valid]))
    time_minutes = (times - t_ref) * 1440.0
    span = max(float(np.ptp(time_minutes[valid])), 1.0)
    cadence = max(float(np.nanmedian(np.diff(np.sort(time_minutes[valid])))), 0.1)
    t0 = (float(initial_midpoint) - t_ref) * 1440.0 if initial_midpoint is not None else 0.0
    duration = float(initial_duration_minutes or max(3.0 * cadence, span * 0.2))
    depth = float(np.clip(initial_depth or 0.01, 1e-5, 0.25))
    ror = float(np.clip(np.sqrt(depth), 0.005, 0.5))
    median_error = max(float(np.nanmedian(error[valid])), _robust_sigma(flux[valid] - np.nanmedian(flux[valid])), 1e-5)
    fit_error = np.maximum(error, median_error * 0.35)

    def model(parameters):
        center, log_duration, log_ror, impact, baseline, slope = parameters
        planet_radius = float(np.exp(log_ror))
        transit_duration = float(np.exp(log_duration))
        contact = max((1.0 + planet_radius) ** 2 - impact**2, 1e-8)
        velocity = 2.0 * np.sqrt(contact) / max(transit_duration, cadence * 0.25)
        separation = np.sqrt(impact**2 + np.square(velocity * (time_minutes - center)))
        delta = quad_limbdark_light_curve(
            float(limb_darkening[0]),
            float(limb_darkening[1]),
            np.ascontiguousarray(separation, dtype=np.float64),
            planet_radius,
        )
        baseline_model = baseline + slope * (time_minutes - center) / span
        return baseline_model + np.asarray(delta, dtype=np.float64)

    def residual(parameters):
        return (model(parameters)[valid] - flux[valid]) / fit_error[valid]

    initial = np.asarray([t0, np.log(duration), np.log(ror), 0.3, 1.0, 0.0])
    lower = np.asarray([
        float(np.min(time_minutes[valid])),
        np.log(max(cadence * 0.5, 0.1)),
        np.log(0.003),
        0.0,
        0.8,
        -0.2,
    ])
    upper = np.asarray([
        float(np.max(time_minutes[valid])),
        np.log(max(span * 0.8, cadence)),
        np.log(0.6),
        1.2,
        1.2,
        0.2,
    ])
    initial = np.clip(initial, lower + 1e-8, upper - 1e-8)
    try:
        solution = least_squares(residual, initial, bounds=(lower, upper), loss="soft_l1", max_nfev=800)
        fitted = model(solution.x)
        chi2 = float(np.sum(np.square(residual(solution.x))))
        null = np.full(np.count_nonzero(valid), np.average(flux[valid], weights=1.0 / np.square(fit_error[valid])))
        null_chi2 = float(np.sum(np.square((flux[valid] - null) / fit_error[valid])))
        n = int(np.count_nonzero(valid))
        delta_bic = float((null_chi2 + np.log(n)) - (chi2 + len(solution.x) * np.log(n)))
        center, log_duration, log_ror, impact, _baseline, _slope = solution.x
        fitted_depth = float(max(0.0, 1.0 - np.min(fitted)))
        snr = fitted_depth / max(median_error / np.sqrt(max(2, n // 4)), 1e-9)
        return TransitModelFit(
            success=bool(solution.success and np.all(np.isfinite(fitted))),
            model_flux=[float(value) for value in fitted],
            mid_transit_jd=float(t_ref + center / 1440.0),
            duration_minutes=float(np.exp(log_duration)),
            radius_ratio=float(np.exp(log_ror)),
            impact_parameter=float(impact),
            depth=fitted_depth,
            snr=float(snr),
            reduced_chi2=float(chi2 / max(n - len(solution.x), 1)),
            delta_bic=delta_bic,
            message="Quadratic limb-darkened transit fit",
        )
    except Exception as exc:
        return TransitModelFit(False, empty_model.tolist(), message=str(exc))


def _weighted_comparison_ensemble(comparisons, uncertainties):
    values = np.asarray(comparisons, dtype=np.float64)
    errors = np.asarray(uncertainties, dtype=np.float64)
    normalized = np.full_like(values, np.nan)
    normalized_errors = np.full_like(values, np.nan)
    stability = np.full(values.shape[1], np.inf, dtype=np.float64)
    flags = []
    for column in range(values.shape[1]):
        median = _positive_median(values[:, column])
        normalized[:, column] = values[:, column] / median
        normalized_errors[:, column] = errors[:, column] / median
        finite = np.isfinite(normalized[:, column]) & (normalized[:, column] > 0)
        if int(np.count_nonzero(finite)) >= 4:
            stability[column] = max(_robust_sigma(normalized[finite, column] - 1.0), 1e-5)
    raw_weights = 1.0 / np.square(stability)
    raw_weights[~np.isfinite(raw_weights)] = 0.0
    if np.sum(raw_weights) <= 0:
        raise ValueError("Karşılaştırma yıldızlarından kararlı bir ensemble üretilemedi.")
    median_weight = float(np.median(raw_weights[raw_weights > 0]))
    raw_weights = np.minimum(raw_weights, median_weight * 25.0)
    weights = raw_weights / np.sum(raw_weights)
    for index, weight in enumerate(weights):
        if weight < 0.02 and values.shape[1] > 1:
            flags.append(f"comparison_{index + 1}_unstable")
    ensemble = np.full(values.shape[0], np.nan, dtype=np.float64)
    ensemble_error = np.full(values.shape[0], np.nan, dtype=np.float64)
    for row in range(values.shape[0]):
        finite = np.isfinite(normalized[row]) & (normalized[row] > 0) & (weights > 0)
        if not np.any(finite):
            continue
        local_weights = weights[finite] / np.sum(weights[finite])
        ensemble[row] = float(np.sum(local_weights * normalized[row, finite]))
        measurement_error = float(np.sqrt(np.sum(np.square(local_weights * normalized_errors[row, finite]))))
        ensemble_error[row] = max(measurement_error, float(np.sqrt(np.sum(np.square(local_weights * stability[finite])))))
    return ensemble, ensemble_error, weights, flags


def _robust_polynomial_trend(times, flux, uncertainty, valid, order, *, exclude=None):
    centered = (times - np.nanmedian(times[valid])) * 1440.0
    scale = max(float(np.ptp(centered[valid])) / 2.0, 1.0)
    x = centered / scale
    excluded = np.asarray(exclude if exclude is not None else np.zeros(len(flux), dtype=bool), dtype=bool)
    keep = valid & ~excluded
    if int(np.count_nonzero(keep)) < 4:
        keep = valid.copy()
    requested_degree = 0 if int(np.count_nonzero(valid)) < 10 else int(order)
    degree = max(0, min(requested_degree, 2, int(np.count_nonzero(keep)) - 2))
    trend = np.full(len(flux), np.nanmedian(flux[valid]), dtype=np.float64)
    for _ in range(4):
        if int(np.count_nonzero(keep)) <= degree + 1:
            break
        weights = 1.0 / np.maximum(uncertainty[keep], np.nanmedian(uncertainty[valid]) * 0.25)
        coefficients = np.polyfit(x[keep], flux[keep], degree, w=weights)
        trend = np.polyval(coefficients, x)
        residual = flux - trend
        sigma = _robust_sigma(residual[keep])
        keep = valid & ~excluded & (residual > -2.5 * sigma) & (residual < 4.0 * sigma)
    fallback = float(np.nanmedian(flux[valid]))
    return np.where(np.isfinite(trend) & (trend > 0), trend, fallback)


def _search_single_transit(times, flux, uncertainty, valid):
    indices = np.flatnonzero(valid)
    ordered = indices[np.argsort(times[indices])]
    values = flux[ordered]
    errors = uncertainty[ordered]
    n = len(ordered)
    best = None
    max_width = max(2, min(n // 2, max(3, n // 3)))
    for width in range(2, max_width + 1):
        for start in range(1, n - width):
            end = start + width
            inside = np.arange(start, end)
            outside = np.concatenate((np.arange(0, start), np.arange(end, n)))
            if len(outside) < 3:
                continue
            baseline = float(np.nanmedian(values[outside]))
            in_level = float(np.nanmedian(values[inside]))
            depth = max(0.0, baseline - in_level)
            noise = max(_robust_sigma(values[outside] - baseline), float(np.nanmedian(errors[outside])), 1e-6)
            depth_uncertainty = float(np.sqrt(noise**2 / len(inside) + noise**2 / len(outside)))
            snr = depth / max(depth_uncertainty, 1e-9)
            coverage_ok = start >= 2 and n - end >= 2
            score = snr * (1.0 if coverage_ok else 0.65)
            baseline_balance = min(start, n - end)
            better_score = best is None or score > best["score"] + 1e-9
            better_tie = (
                best is not None
                and abs(score - best["score"]) <= 1e-9
                and baseline_balance > best["baseline_balance"]
            )
            if better_score or better_tie:
                cadence = float(np.nanmedian(np.diff(times[ordered]))) * 1440.0 if n > 1 else 0.0
                duration = (times[ordered[end - 1]] - times[ordered[start]]) * 1440.0 + max(cadence, 0.0)
                best = {
                    "score": float(score),
                    "baseline_balance": int(baseline_balance),
                    "depth": float(depth),
                    "depth_uncertainty": depth_uncertainty,
                    "snr": float(snr),
                    "midpoint_jd": float(np.mean(times[ordered[inside]])),
                    "duration_minutes": float(max(duration, cadence)),
                    "in_transit_points": int(width),
                    "in_transit_indices": [int(value) for value in ordered[inside]],
                    "coverage_ok": bool(coverage_ok),
                }
    if best is None:
        return {
            "score": 0.0,
            "baseline_balance": 0,
            "depth": 0.0,
            "depth_uncertainty": float("inf"),
            "snr": 0.0,
            "midpoint_jd": float(np.nanmedian(times[valid])),
            "duration_minutes": 0.0,
            "in_transit_points": 0,
            "in_transit_indices": [],
            "coverage_ok": False,
        }
    return best


def _coerce_uncertainties(values, flux):
    if values is not None:
        array = np.asarray(values, dtype=np.float64)
        if array.shape == flux.shape:
            finite = array[np.isfinite(array) & (array > 0)]
            fallback = float(np.nanmedian(finite)) if finite.size else max(_robust_sigma(flux), 1.0)
            return np.where(np.isfinite(array) & (array > 0), array, fallback)
    scatter = max(_robust_sigma(flux - np.nanmedian(flux)), _positive_median(flux) * 0.001, 1.0)
    return np.full_like(flux, scatter, dtype=np.float64)


def _coerce_comparison_uncertainties(values, comparisons):
    if values is not None:
        array = np.asarray(values, dtype=np.float64)
        if array.ndim == 1:
            array = array[:, None]
        if array.shape == comparisons.shape:
            result = array.copy()
            for column in range(result.shape[1]):
                finite = result[:, column][np.isfinite(result[:, column]) & (result[:, column] > 0)]
                fallback = float(np.nanmedian(finite)) if finite.size else max(_robust_sigma(comparisons[:, column]), 1.0)
                result[:, column] = np.where(np.isfinite(result[:, column]) & (result[:, column] > 0), result[:, column], fallback)
            return result
    result = np.empty_like(comparisons)
    for column in range(comparisons.shape[1]):
        scatter = max(_robust_sigma(comparisons[:, column]), _positive_median(comparisons[:, column]) * 0.001, 1.0)
        result[:, column] = scatter
    return result


def _positive_median(values):
    array = np.asarray(values, dtype=np.float64)
    valid = array[np.isfinite(array) & (array > 0)]
    if valid.size == 0:
        raise ValueError("Pozitif fotometri ölçümü bulunamadı.")
    return float(np.nanmedian(valid))


def _robust_sigma(values):
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return 1e-6
    median = float(np.nanmedian(finite))
    mad = float(np.nanmedian(np.abs(finite - median)))
    return max(1.4826 * mad, float(np.nanstd(finite)) * 0.15, 1e-6)
