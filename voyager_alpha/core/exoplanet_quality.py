from __future__ import annotations

import json
from dataclasses import dataclass, field
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

from .calibration import calibrate_science_frame, load_master_dark, load_master_frame
from .exoplanet import ApertureMeasurement, aperture_measurement
from .fits_io import read_fits_header, read_fits_image
from .metadata import inspect_sequence
from .registration import detect_registration_stars, register_frame
from .models import RegistrationSolution


LINEARITY_KEYS = ("LINLEVEL", "LINEAR", "LINEARIT", "LINEARITY")
SATURATION_KEYS = ("SATURATE", "SATLEVEL", "SATUR", "MAXADU")


@dataclass
class DetectorLinearity:
    limit_adu: float | None
    source: str
    verified: bool
    digital_max_adu: float | None = None


@dataclass
class StarQualityAssessment:
    x: float
    y: float
    role: str
    status: str
    median_flux: float
    median_snr: float
    peak_adu: float
    peak_fraction: float | None
    median_fwhm_px: float
    valid_fraction: float
    stability: float | None = None
    delta_mag: float | None = None
    nearest_neighbor_px: float | None = None
    gaia_g_mag: float | None = None
    bp_rp: float | None = None
    color_delta: float | None = None
    variable_status: str = "unchecked"
    flags: list[str] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)


@dataclass
class PhotometryPreflightResult:
    linearity: DetectorLinearity
    target: StarQualityAssessment
    references: list[StarQualityAssessment]
    recommended_xy: list[tuple[float, float]]
    warnings: list[str] = field(default_factory=list)
    aligned_frames: int = 0
    total_frames: int = 0

    @property
    def analysis_allowed(self) -> bool:
        return self.target.status != "FAIL" and len(self.recommended_xy) >= 2


def resolve_detector_linearity(header, override_adu: float = 0.0) -> DetectorLinearity:
    if np.isfinite(override_adu) and float(override_adu) > 0:
        return DetectorLinearity(float(override_adu), "user", True, _digital_max_from_header(header))
    for key in LINEARITY_KEYS:
        value = _positive_header_float(header, key)
        if value is not None:
            return DetectorLinearity(value, f"FITS:{key}", True, _digital_max_from_header(header))
    for key in SATURATION_KEYS:
        value = _positive_header_float(header, key)
        if value is not None:
            return DetectorLinearity(value * 0.95, f"FITS:{key} x 0.95", False, value)
    digital_max = _digital_max_from_header(header)
    if digital_max is not None:
        return DetectorLinearity(digital_max * 0.90, "FITS numeric range x 0.90", False, digital_max)
    return DetectorLinearity(None, "unknown", False, None)


def run_photometry_preflight(
    file_paths,
    target_xy,
    *,
    aperture_radius: float = 6.0,
    linearity_limit_adu: float = 0.0,
    master_paths=None,
    reference_header=None,
    catalog_lookup: bool = True,
    max_references: int = 8,
    progress_callback=None,
    stop_callback=None,
) -> PhotometryPreflightResult:
    records = inspect_sequence(list(file_paths))
    if len(records) < 5:
        raise ValueError("Fotometri ön kontrolü için en az 5 FITS karesi gerekir.")
    records.sort(key=lambda record: record.midpoint_jd or 0.0)
    header = read_fits_header(records[0].file_path)
    linearity = resolve_detector_linearity(header, linearity_limit_adu)
    if linearity.limit_adu is None:
        raise ValueError(
            "Kamera lineerlik sınırı FITS başlığından belirlenemedi. "
            "Kamera lineerlik testinde ölçülen ADU sınırını Photometry Settings alanına girin."
        )

    masters = _load_masters(master_paths or {})
    reference_raw = np.asarray(read_fits_image(records[0].file_path), dtype=np.float32)
    reference = _calibrated_frame(reference_raw, records[0], masters)
    annulus_outer = max(aperture_radius * 2.3, aperture_radius + 7.0)
    detected = detect_registration_stars(
        reference,
        threshold_sigma=6.0,
        max_stars=180,
        edge_margin=max(12, int(np.ceil(annulus_outer + 3.0))),
    )
    target_x, target_y = map(float, target_xy)
    target_reference = aperture_measurement(
        reference,
        target_x,
        target_y,
        aperture_radius=aperture_radius,
        gain_e_per_adu=_usable_gain(records[0]),
    )
    target_reference = merge_detector_quality(
        target_reference,
        aperture_measurement(
            reference_raw,
            target_x,
            target_y,
            aperture_radius=aperture_radius,
            linearity_limit_adu=linearity.limit_adu,
        ),
    )
    candidate_points = _initial_candidate_points(
        reference,
        reference_raw,
        detected,
        target_reference,
        aperture_radius,
        annulus_outer,
        linearity.limit_adu,
        _usable_gain(records[0]),
    )
    if len(candidate_points) < 2:
        raise ValueError("Hedef çevresinde fotometriye uygun en az iki referans yıldızı bulunamadı.")

    target_measurements: list[ApertureMeasurement] = []
    candidate_measurements: list[list[ApertureMeasurement]] = [[] for _ in candidate_points]
    aligned_count = 0
    total = len(records)
    for index, record in enumerate(records):
        if stop_callback and stop_callback():
            raise RuntimeError("Fotometri ön kontrolü durduruldu.")
        raw = np.asarray(read_fits_image(record.file_path), dtype=np.float32)
        data = _calibrated_frame(raw, record, masters)
        if index == 0:
            aligned = reference
            solution = RegistrationSolution.identity()
            reliable = True
        else:
            aligned, solution = register_frame(reference, data)
            reliable = solution.method.startswith("star") and solution.rms_px <= 2.0
        if reliable:
            aligned_count += 1
            target_measurements.append(
                _science_and_detector_measurement(
                    aligned,
                    raw,
                    (target_x, target_y),
                    solution,
                    aperture_radius,
                    linearity.limit_adu,
                    _usable_gain(record),
                )
            )
            for column, point in enumerate(candidate_points):
                candidate_measurements[column].append(
                    _science_and_detector_measurement(
                        aligned,
                        raw,
                        point,
                        solution,
                        aperture_radius,
                        linearity.limit_adu,
                        _usable_gain(record),
                    )
                )
        else:
            target_measurements.append(_invalid_measurement(target_x, target_y, "weak_registration"))
            for column, point in enumerate(candidate_points):
                candidate_measurements[column].append(_invalid_measurement(*point, "weak_registration"))
        if progress_callback:
            progress_callback(10 + int(65 * (index + 1) / total), f"Yıldız kalitesi {index + 1}/{total}")

    catalog = _catalog_annotations(
        reference_header or header,
        [(target_x, target_y), *candidate_points],
        enabled=catalog_lookup,
    )
    target = _assess_star(
        target_x,
        target_y,
        "target",
        target_measurements,
        linearity,
        nearest_neighbor_px=_nearest_neighbor((target_x, target_y), detected),
        catalog=catalog.get(0),
    )
    target_flux = max(target.median_flux, 1e-12)
    raw_references = []
    for index, (point, measurements) in enumerate(zip(candidate_points, candidate_measurements), start=1):
        median_flux = _positive_median([item.flux for item in measurements])
        delta_mag = float(-2.5 * np.log10(max(median_flux, 1e-12) / target_flux))
        assessment = _assess_star(
            *point,
            "reference",
            measurements,
            linearity,
            delta_mag=delta_mag,
            nearest_neighbor_px=_nearest_neighbor(point, detected),
            catalog=catalog.get(index),
        )
        raw_references.append(assessment)

    _apply_leave_one_out_stability(raw_references, candidate_measurements)
    target_color = target.bp_rp
    for assessment in raw_references:
        if target_color is not None and assessment.bp_rp is not None:
            assessment.color_delta = abs(assessment.bp_rp - target_color)
            if assessment.color_delta > 0.7:
                assessment.flags.append("large_color_difference")
                assessment.messages.append(f"Gaia BP-RP renk farkı {assessment.color_delta:.2f} mag")
                if assessment.status == "PASS":
                    assessment.status = "CAUTION"

    recommended = _select_reference_ensemble(raw_references, max_references=max_references)
    warnings = []
    if not linearity.verified:
        warnings.append(
            "Lineerlik sınırı kamera testiyle doğrulanmadı; FITS sayısal aralığından konservatif olarak türetildi."
        )
    if aligned_count < total:
        warnings.append(f"{total - aligned_count} kare güvenilir yıldız hizası kuramadı.")
    if catalog_lookup and not any(item.variable_status != "unchecked" for item in raw_references):
        warnings.append("Gaia/VSX değişkenlik kontrolü tamamlanamadı; sekans içi kararlılık kullanıldı.")
    if len(recommended) < 3:
        warnings.append("Yalnız iki uygun referans bulundu; bilimsel hedef en az üç, ideal olarak 8 referanstır.")
    if progress_callback:
        progress_callback(100, f"Ön kontrol tamamlandı: {len(recommended)} referans")
    return PhotometryPreflightResult(
        linearity=linearity,
        target=target,
        references=raw_references,
        recommended_xy=[(item.x, item.y) for item in recommended],
        warnings=warnings,
        aligned_frames=aligned_count,
        total_frames=total,
    )


def _load_masters(master_paths):
    dark = None
    dark_exposure = None
    if master_paths.get("dark"):
        dark, dark_exposure = load_master_dark(master_paths["dark"])
    return {
        "bias": load_master_frame(master_paths["bias"]) if master_paths.get("bias") else None,
        "dark": dark,
        "dark_exposure": dark_exposure,
        "flat": load_master_frame(master_paths["flat"]) if master_paths.get("flat") else None,
    }


def _calibrated_frame(raw, record, masters):
    return calibrate_science_frame(
        raw,
        master_bias=masters["bias"],
        master_dark=masters["dark"],
        master_flat=masters["flat"],
        science_exposure=record.exposure_seconds,
        dark_exposure=masters["dark_exposure"],
    )


def _initial_candidate_points(
    image,
    raw_image,
    detected,
    target,
    aperture_radius,
    annulus_outer,
    linearity_limit,
    gain,
):
    candidates = []
    target_flux = max(float(target.flux), 1e-12)
    for x, y in detected:
        if np.hypot(x - target.x, y - target.y) <= annulus_outer * 2.0:
            continue
        measurement = aperture_measurement(
            image,
            float(x),
            float(y),
            aperture_radius=aperture_radius,
            gain_e_per_adu=gain,
        )
        measurement = merge_detector_quality(
            measurement,
            aperture_measurement(
                raw_image,
                float(x),
                float(y),
                aperture_radius=aperture_radius,
                linearity_limit_adu=linearity_limit,
            ),
        )
        if not np.isfinite(measurement.flux) or measurement.flux <= 0:
            continue
        if "saturated" in measurement.flags or "outside_frame" in measurement.flags:
            continue
        delta_mag = -2.5 * np.log10(measurement.flux / target_flux)
        nearest = _nearest_neighbor((measurement.x, measurement.y), detected)
        if nearest is not None and nearest < max(annulus_outer, measurement.fwhm_px * 2.5):
            continue
        if -1.5 <= delta_mag <= 1.5:
            candidates.append((abs(delta_mag), -measurement.snr, float(measurement.x), float(measurement.y)))
    candidates.sort()
    return [(x, y) for _delta, _snr, x, y in candidates[:36]]


def _assess_star(
    x,
    y,
    role,
    measurements,
    linearity,
    *,
    delta_mag=None,
    nearest_neighbor_px=None,
    catalog=None,
):
    valid = [item for item in measurements if np.isfinite(item.flux) and item.flux > 0]
    valid_fraction = len(valid) / max(len(measurements), 1)
    flags = sorted({flag for item in measurements for flag in item.flags})
    messages = []
    status = "PASS"
    saturated_count = sum("saturated" in item.flags for item in measurements)
    near_count = sum("near_linearity_limit" in item.flags for item in measurements)
    if saturated_count:
        status = "FAIL"
        messages.append(f"{saturated_count} karede lineerlik sınırı aşıldı")
    if valid_fraction < 0.8:
        status = "FAIL"
        messages.append(f"Geçerli ölçüm oranı %{valid_fraction * 100:.0f}")
    elif valid_fraction < 0.95 and status != "FAIL":
        status = "CAUTION"
        messages.append(f"Geçerli ölçüm oranı %{valid_fraction * 100:.0f}")
    if near_count and status == "PASS":
        status = "CAUTION"
        messages.append(f"{near_count} kare lineerlik sınırının %90 üzerinde")
    if not linearity.verified and status == "PASS":
        status = "CAUTION"
        flags.append("linearity_limit_inferred")
        messages.append("Kamera lineerlik sınırı ölçüm profiliyle doğrulanmadı")
    median_fwhm = _median([item.fwhm_px for item in valid])
    median_snr = _median([item.snr for item in valid])
    peak = max((item.peak_adu for item in valid), default=float("nan"))
    peak_fraction = peak / linearity.limit_adu if linearity.limit_adu and np.isfinite(peak) else None
    if np.isfinite(median_fwhm) and median_fwhm < 2.0:
        flags.append("undersampled_psf")
        messages.append(f"PSF FWHM {median_fwhm:.2f} px; hafif defocus önerilir")
        if status == "PASS":
            status = "CAUTION"
    elif np.isfinite(median_fwhm) and median_fwhm > 12.0:
        flags.append("broad_psf")
        messages.append(f"PSF FWHM {median_fwhm:.1f} px; crowding ve aperture kontrol edilmeli")
        if status == "PASS":
            status = "CAUTION"
    if np.isfinite(median_snr) and median_snr < 50:
        flags.append("low_snr")
        messages.append(f"Medyan SNR {median_snr:.1f}; hedef en az 50")
        if status == "PASS":
            status = "CAUTION"
    crowding_limit = max(10.0, median_fwhm * 2.5) if np.isfinite(median_fwhm) else 10.0
    if nearest_neighbor_px is not None and nearest_neighbor_px < crowding_limit:
        flags.append("crowded_aperture")
        messages.append(f"En yakın yıldız {nearest_neighbor_px:.1f} px; aperture karışımı riski")
        status = "FAIL"
    gaia_g = None
    bp_rp = None
    variable_status = "unchecked"
    if catalog:
        gaia_g = catalog.get("g_mag")
        bp_rp = catalog.get("bp_rp")
        variable_status = catalog.get("variable_status", "unchecked")
        if variable_status == "variable":
            flags.append("catalog_variable")
            messages.append("Gaia/VSX değişken yıldız kaydı")
            if role == "reference":
                status = "FAIL"
    return StarQualityAssessment(
        x=float(x),
        y=float(y),
        role=role,
        status=status,
        median_flux=_positive_median([item.flux for item in valid]),
        median_snr=median_snr,
        peak_adu=peak,
        peak_fraction=peak_fraction,
        median_fwhm_px=median_fwhm,
        valid_fraction=valid_fraction,
        delta_mag=delta_mag,
        nearest_neighbor_px=nearest_neighbor_px,
        gaia_g_mag=gaia_g,
        bp_rp=bp_rp,
        variable_status=variable_status,
        flags=sorted(set(flags)),
        messages=messages,
    )


def _apply_leave_one_out_stability(assessments, measurement_columns):
    if not assessments:
        return
    values = np.asarray(
        [[measurement.flux for measurement in column] for column in measurement_columns],
        dtype=np.float64,
    ).T
    medians = np.asarray([_positive_median(values[:, column]) for column in range(values.shape[1])])
    normalized = values / np.maximum(medians, 1e-12)
    for column, assessment in enumerate(assessments):
        others = np.delete(normalized, column, axis=1)
        if others.shape[1] == 0:
            assessment.stability = float("inf")
            assessment.status = "FAIL"
            assessment.flags.append("no_ensemble_check")
            continue
        ensemble = np.nanmedian(others, axis=1)
        ratio = normalized[:, column] / ensemble
        finite = np.isfinite(ratio) & (ratio > 0)
        assessment.stability = _robust_sigma(ratio[finite] - 1.0) if np.count_nonzero(finite) >= 4 else float("inf")
    finite_stability = np.asarray(
        [item.stability for item in assessments if item.stability is not None and np.isfinite(item.stability)]
    )
    baseline = float(np.median(finite_stability)) if finite_stability.size else 0.005
    limit = max(0.005, baseline * 3.0)
    for assessment in assessments:
        if assessment.stability is None or not np.isfinite(assessment.stability) or assessment.stability > limit:
            assessment.flags.append("unstable_reference")
            assessment.messages.append(
                "Sekans içi kararlılık yetersiz"
                if not np.isfinite(assessment.stability or np.inf)
                else f"Leave-one-out saçılım %{assessment.stability * 100:.3f}"
            )
            assessment.status = "FAIL"


def _select_reference_ensemble(assessments, *, max_references):
    usable = [item for item in assessments if item.status != "FAIL"]
    ideal = [item for item in usable if item.delta_mag is not None and -0.44 <= item.delta_mag <= 0.75]
    pool = ideal if len(ideal) >= 2 else usable
    pool.sort(
        key=lambda item: (
            item.status != "PASS",
            item.variable_status == "unchecked",
            item.color_delta if item.color_delta is not None else 0.5,
            item.stability if item.stability is not None else float("inf"),
            abs(item.delta_mag or 0.0),
        )
    )
    return pool[: max(2, int(max_references))]


def _catalog_annotations(header, points, *, enabled):
    if not enabled:
        return {}
    try:
        wcs = WCS(header).celestial
        if not wcs.has_celestial:
            return {}
        pixels = np.asarray(points, dtype=np.float64)
        world = wcs.all_pix2world(pixels, 0)
        coords = SkyCoord(world[:, 0] * u.deg, world[:, 1] * u.deg)
    except Exception:
        return {}
    result = {index: {"variable_status": "unchecked"} for index in range(len(points))}
    _apply_gaia_annotations(coords, result)
    _apply_vsx_annotations(coords, result)
    return result


def _apply_gaia_annotations(coords, result):
    try:
        from astroquery.gaia import Gaia

        center = SkyCoord(np.mean(coords.ra.deg) * u.deg, np.mean(coords.dec.deg) * u.deg)
        radius = max(float(np.max(center.separation(coords).deg)) + 0.01, 0.02)
        Gaia.TIMEOUT = 15
        query = (
            "SELECT TOP 500 source_id,ra,dec,phot_g_mean_mag,bp_rp,phot_variable_flag "
            "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
            f"CIRCLE('ICRS',{center.ra.deg:.9f},{center.dec.deg:.9f},{radius:.7f}))"
        )
        table = Gaia.launch_job_async(query, dump_to_file=False).get_results()
        if len(table) == 0:
            return
        catalog_coords = SkyCoord(np.asarray(table["ra"], dtype=float) * u.deg, np.asarray(table["dec"], dtype=float) * u.deg)
        nearest, separation, _ = coords.match_to_catalog_sky(catalog_coords)
        for index, (row, sep) in enumerate(zip(nearest, separation.arcsec)):
            if sep > 3.0:
                continue
            g_mag = _finite_or_none(table["phot_g_mean_mag"][row])
            bp_rp = _finite_or_none(table["bp_rp"][row])
            variable = str(table["phot_variable_flag"][row]).strip().upper()
            result[index].update(
                g_mag=g_mag,
                bp_rp=bp_rp,
                variable_status="variable" if variable == "VARIABLE" else "not_flagged",
            )
    except Exception:
        return


def _apply_vsx_annotations(coords, result):
    try:
        center = SkyCoord(np.mean(coords.ra.deg) * u.deg, np.mean(coords.dec.deg) * u.deg)
        radius = max(float(np.max(center.separation(coords).deg)) + 0.01, 0.02)
        query = urlencode(
            {
                "view": "api.list",
                "ra": f"{center.ra.deg:.8f}",
                "dec": f"{center.dec.deg:.8f}",
                "radius": f"{radius:.6f}",
                "format": "json",
            }
        )
        with urlopen(f"https://www.aavso.org/vsx/index.php?{query}", timeout=12) as response:
            payload = json.load(response)
        rows = payload.get("VSXObjects", {}).get("VSXObject", [])
        if isinstance(rows, dict):
            rows = [rows]
        vsx_coords = []
        for row in rows:
            try:
                vsx_coords.append(
                    SkyCoord(
                        f"{row['RA2000']} {row['Declination2000']}",
                        unit=(u.hourangle, u.deg),
                    )
                )
            except Exception:
                continue
        if not vsx_coords:
            return
        catalog = SkyCoord(vsx_coords)
        nearest, separation, _ = coords.match_to_catalog_sky(catalog)
        for index, sep in enumerate(separation.arcsec):
            if sep <= 5.0:
                result[index]["variable_status"] = "variable"
    except Exception:
        return


def _nearest_neighbor(point, stars):
    if len(stars) < 2:
        return None
    distances = np.hypot(stars[:, 0] - point[0], stars[:, 1] - point[1])
    positive = distances[distances > 2.0]
    return float(np.min(positive)) if positive.size else None


def _invalid_measurement(x, y, flag):
    return ApertureMeasurement(
        float("nan"),
        float("nan"),
        float(x),
        float(y),
        0.0,
        0.0,
        0,
        [flag],
    )


def sensor_coordinates(point, solution: RegistrationSolution) -> tuple[float, float]:
    aligned = np.asarray(point, dtype=np.float64)
    try:
        sensor = np.linalg.solve(solution.matrix_xy, aligned - solution.offset_xy)
        if np.all(np.isfinite(sensor)):
            return float(sensor[0]), float(sensor[1])
    except (ValueError, np.linalg.LinAlgError):
        pass
    return float(point[0]), float(point[1])


def merge_detector_quality(
    science: ApertureMeasurement,
    detector: ApertureMeasurement,
) -> ApertureMeasurement:
    science.peak_adu = detector.peak_adu
    science.linearity_limit_adu = detector.linearity_limit_adu
    science.saturated_pixels = detector.saturated_pixels
    for flag in ("saturated", "near_linearity_limit"):
        if flag in detector.flags and flag not in science.flags:
            science.flags.append(flag)
    return science


def _science_and_detector_measurement(
    aligned,
    raw,
    point,
    solution,
    aperture_radius,
    linearity_limit,
    gain,
):
    science = aperture_measurement(
        aligned,
        *point,
        aperture_radius=aperture_radius,
        gain_e_per_adu=gain,
    )
    sensor_xy = sensor_coordinates(point, solution)
    detector = aperture_measurement(
        raw,
        *sensor_xy,
        aperture_radius=aperture_radius,
        linearity_limit_adu=linearity_limit,
    )
    return merge_detector_quality(science, detector)


def _usable_gain(record):
    gain = record.camera.gain
    if gain is None or gain <= 0:
        return None
    if record.camera.gain_keyword == "EGAIN" or gain <= 10.0:
        return float(gain)
    return None


def _digital_max_from_header(header):
    try:
        bitpix = int(header.get("BITPIX", 0))
        if bitpix <= 0 or bitpix > 32:
            return None
        bscale = float(header.get("BSCALE", 1.0))
        bzero = float(header.get("BZERO", 0.0))
        return float((2 ** (bitpix - 1) - 1) * bscale + bzero)
    except (TypeError, ValueError, OverflowError):
        return None


def _positive_header_float(header, key):
    try:
        value = float(header.get(key))
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


def _positive_median(values):
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array) & (array > 0)]
    return float(np.median(finite)) if finite.size else float("nan")


def _median(values):
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    return float(np.median(finite)) if finite.size else float("nan")


def _robust_sigma(values):
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return float("inf")
    median = float(np.median(finite))
    return float(max(1.4826 * np.median(np.abs(finite - median)), 1e-8))


def _finite_or_none(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None
