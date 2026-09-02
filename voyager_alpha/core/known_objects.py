from __future__ import annotations

import warnings
import csv
import io
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS

from .detection import robust_location_scale
from .discovery_method import DOCUMENTED_DISCOVERY_METHOD
from .astrometry import estimate_pixel_scale_arcsec
from .models import Detection, FrameRecord, KnownObjectPrediction, Tracklet


@dataclass
class KnownObject:
    name: str
    ra: float
    dec: float
    magnitude: float | None
    raw: dict[str, Any]
    motion_ra_arcsec_hour: float | None = None
    motion_dec_arcsec_hour: float | None = None
    object_type: str = "asteroid"


class KnownObjectMatcher:
    def __init__(self, tolerance_arcsec: float = 20.0):
        self.tolerance_arcsec = float(tolerance_arcsec)
        self._cache: dict[tuple[float, float, float, int], list[KnownObject]] = {}

    def objects_for_frame(self, frame: FrameRecord, header) -> list[KnownObject]:
        center, radius = field_center_and_radius(header, frame.shape)
        if center is None or radius is None or frame.midpoint_jd is None:
            return []
        cache_key = (
            round(float(center.ra.degree), 3),
            round(float(center.dec.degree), 3),
            round(float(radius), 3),
            int(round(float(frame.midpoint_jd) * 1440.0)),
        )
        if cache_key not in self._cache:
            self._cache[cache_key] = self._query(center, radius, frame.midpoint_jd)
        return self._cache[cache_key]

    def predictions_for_frame(
        self,
        frame: FrameRecord,
        header,
        image: np.ndarray | None = None,
    ) -> list[KnownObjectPrediction]:
        objects = self.objects_for_frame(frame, header)
        wcs = WCS(header)
        pixel_scale = estimate_pixel_scale_arcsec(header)
        predictions = []
        for item in objects:
            prediction = KnownObjectPrediction(
                name=item.name,
                ra=item.ra,
                dec=item.dec,
                magnitude=item.magnitude,
                object_type=item.object_type,
                motion_ra_arcsec_hour=item.motion_ra_arcsec_hour,
                motion_dec_arcsec_hour=item.motion_dec_arcsec_hour,
                raw=item.raw,
            )
            try:
                x, y = wcs.all_world2pix(item.ra, item.dec, 0)
                prediction.x, prediction.y = float(x), float(y)
            except Exception:
                predictions.append(prediction)
                continue
            if image is not None and _inside_image(image, prediction.x, prediction.y):
                motion = float(np.hypot(item.motion_ra_arcsec_hour or 0.0, item.motion_dec_arcsec_hour or 0.0))
                exposure = float(frame.exposure_seconds or 0.0)
                trail_px = motion * exposure / 3600.0 / max(float(pixel_scale or 0.0), 1e-6)
                search_radius = int(np.clip(round(5.0 + 0.5 * trail_px), 5, 18))
                snr, offset = measure_local_peak(
                    image,
                    prediction.x,
                    prediction.y,
                    radius_px=search_radius,
                )
                prediction.local_snr = snr
                prediction.local_offset_px = offset
                prediction.expected_trail_px = float(trail_px)
                prediction.search_radius_px = search_radius
                height, width = image.shape
                prediction.near_edge = bool(
                    prediction.x < search_radius
                    or prediction.y < search_radius
                    or prediction.x >= width - search_radius
                    or prediction.y >= height - search_radius
                )
                prediction.visible = bool(
                    snr >= DOCUMENTED_DISCOVERY_METHOD.known_visible_snr
                    and offset <= DOCUMENTED_DISCOVERY_METHOD.known_visible_offset_px
                )
                prediction.confidence = prediction_confidence(prediction)
                prediction.status = prediction_status(prediction)
            predictions.append(prediction)
        return predictions

    def match_detection(self, detection: Detection, objects: list[KnownObject]) -> dict[str, Any] | None:
        if detection.ra is None or detection.dec is None or not objects:
            return None
        target = SkyCoord(ra=detection.ra * u.degree, dec=detection.dec * u.degree, frame="icrs")
        best = None
        best_sep = self.tolerance_arcsec
        for known in objects:
            known_coord = SkyCoord(ra=known.ra * u.degree, dec=known.dec * u.degree, frame="icrs")
            separation = float(target.separation(known_coord).arcsec)
            if separation <= best_sep:
                best_sep = separation
                best = {
                    "name": known.name,
                    "mag": known.magnitude,
                    "sep": separation,
                    "ra": known.ra,
                    "dec": known.dec,
                    "motion_ra": known.motion_ra_arcsec_hour,
                    "motion_dec": known.motion_dec_arcsec_hour,
                    "type": known.object_type,
                }
        return best

    def match_tracklet(
        self,
        tracklet: Tracklet,
        frame: FrameRecord,
        header,
    ) -> dict[str, Any] | None:
        representative = min(
            tracklet.detections,
            key=lambda detection: abs(detection.frame_index - frame.index),
        )
        return self.match_detection(representative, self.objects_for_frame(frame, header))

    def _query(self, center: SkyCoord, radius_deg: float, midpoint_jd: float) -> list[KnownObject]:
        try:
            from astroquery.imcce import Skybot

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rows = Skybot.cone_search(
                    center,
                    rad=float(radius_deg) * u.degree,
                    epoch=Time(float(midpoint_jd), format="jd", scale="utc"),
                )
        except Exception:
            return []
        if rows is None:
            return []
        objects = []
        for row in rows:
            try:
                raw = {name: _plain_value(row[name]) for name in row.colnames}
                objects.append(
                    KnownObject(
                        name=str(_row_value(row, ("Name", "Num"), "Unknown")),
                        ra=float(_row_value(row, ("RA",))),
                        dec=float(_row_value(row, ("DEC", "Dec"))),
                        magnitude=_safe_float(_row_value(row, ("V", "Mv"))),
                        raw=raw,
                        motion_ra_arcsec_hour=_safe_float(
                            _row_value(row, ("dRA(arcsec/h)", "dRA", "RA_rate"))
                        ),
                        motion_dec_arcsec_hour=_safe_float(
                            _row_value(row, ("dDEC(arcsec/h)", "dDEC", "DEC_rate"))
                        ),
                        object_type=str(_row_value(row, ("Type", "Class"), "asteroid")),
                    )
                )
            except Exception:
                continue
        return objects


def field_center_and_radius(header, shape: tuple[int, int]) -> tuple[SkyCoord | None, float | None]:
    try:
        wcs = WCS(header)
        if not wcs.has_celestial:
            return None, None
        height, width = shape
        center_ra, center_dec = wcs.all_pix2world(width / 2.0, height / 2.0, 0)
        center = SkyCoord(float(center_ra) * u.degree, float(center_dec) * u.degree, frame="icrs")
        corner_pixels = [(0, 0), (width - 1, 0), (0, height - 1), (width - 1, height - 1)]
        radius = 0.0
        for x, y in corner_pixels:
            ra, dec = wcs.all_pix2world(x, y, 0)
            corner = SkyCoord(float(ra) * u.degree, float(dec) * u.degree, frame="icrs")
            radius = max(radius, float(center.separation(corner).degree))
        return center, radius * 1.08
    except Exception:
        return None, None


def measure_local_peak(
    image: np.ndarray,
    predicted_x: float,
    predicted_y: float,
    *,
    radius_px: int = 10,
) -> tuple[float, float]:
    data = np.asarray(image, dtype=np.float32)
    x0, x1 = max(0, int(predicted_x) - radius_px), min(data.shape[1], int(predicted_x) + radius_px + 1)
    y0, y1 = max(0, int(predicted_y) - radius_px), min(data.shape[0], int(predicted_y) + radius_px + 1)
    patch = data[y0:y1, x0:x1]
    if patch.size == 0:
        return 0.0, float("inf")
    median, rms = robust_location_scale(patch)
    peak_index = int(np.nanargmax(patch))
    peak_y, peak_x = np.unravel_index(peak_index, patch.shape)
    x = x0 + peak_x
    y = y0 + peak_y
    snr = (float(patch[peak_y, peak_x]) - median) / max(rms, 1e-6)
    return float(snr), float(np.hypot(x - predicted_x, y - predicted_y))


def prediction_confidence(
    prediction: KnownObjectPrediction,
    limiting_magnitude: float = 18.0,
) -> float:
    score = 0.1
    if prediction.magnitude is not None:
        score += min(
            0.25,
            (float(limiting_magnitude) - float(prediction.magnitude)) / max(2.0, float(limiting_magnitude)),
        )
    if prediction.local_snr is not None:
        score += min(0.4, max(0.0, prediction.local_snr) / 20.0)
    if prediction.local_offset_px is not None:
        score += 0.25 * (1.0 - min(float(prediction.local_offset_px), 10.0) / 10.0)
    if float(prediction.expected_trail_px or 0.0) > 1.5:
        score += 0.05
    if prediction.near_edge:
        score -= 0.15
    return round(float(np.clip(score, 0.0, 0.99)), 3)


def prediction_status(prediction: KnownObjectPrediction) -> str:
    if prediction.visible and prediction.confidence >= 0.75:
        return "high_confidence_match"
    if prediction.visible and prediction.confidence >= 0.5:
        return "plausible_visible_match"
    if prediction.magnitude is not None and prediction.magnitude <= 18.0:
        return "predicted_visual_confirmation_weak"
    return "predicted_in_field"


def estimate_visible_limit(predictions: list[KnownObjectPrediction]) -> float | None:
    visible = [item.magnitude for item in predictions if item.visible and item.magnitude is not None]
    return float(max(visible)) if visible else None


def estimate_gaia_visible_limit(
    image: np.ndarray,
    header,
    *,
    catalog_rows: list[dict[str, float]] | None = None,
    start_magnitude: float = 15.0,
    max_magnitude: float = 22.5,
) -> float | None:
    """Estimate field depth from Gaia G stars using the documented bin-recovery rule."""

    rows = catalog_rows if catalog_rows is not None else _query_gaia_field(header, image.shape)
    if not rows:
        return None
    try:
        wcs = WCS(header)
    except Exception:
        return None
    bin_width = DOCUMENTED_DISCOVERY_METHOD.gaia_bin_width_mag
    successful_limit = None
    consecutive_failures = 0
    lower = float(start_magnitude)
    while lower <= float(max_magnitude) + 1e-6:
        upper = lower + bin_width
        in_bin = sorted(
            (row for row in rows if lower <= float(row["phot_g_mean_mag"]) < upper),
            key=lambda row: float(row["phot_g_mean_mag"]),
        )[: DOCUMENTED_DISCOVERY_METHOD.gaia_samples_per_bin]
        recovered = 0
        recovered_magnitudes = []
        for row in in_bin:
            try:
                x, y = wcs.all_world2pix(float(row["ra"]), float(row["dec"]), 0)
            except Exception:
                continue
            if not _inside_image(image, float(x), float(y)):
                continue
            snr, offset = measure_local_peak(image, float(x), float(y), radius_px=5)
            if snr >= 4.0 and offset <= 5.0:
                recovered += 1
                recovered_magnitudes.append(float(row["phot_g_mean_mag"]))
        if recovered >= DOCUMENTED_DISCOVERY_METHOD.gaia_required_recoveries:
            successful_limit = max(recovered_magnitudes)
            consecutive_failures = 0
        else:
            consecutive_failures += 1
            if successful_limit is not None and consecutive_failures >= 2:
                break
        lower = upper
    return float(successful_limit) if successful_limit is not None else None


def _query_gaia_field(header, shape: tuple[int, int], timeout_seconds: int = 12) -> list[dict[str, float]]:
    center, radius = field_center_and_radius(header, shape)
    if center is None or radius is None:
        return []
    query = (
        "SELECT TOP 3000 ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source "
        "WHERE phot_g_mean_mag BETWEEN 14.5 AND 23.0 AND "
        "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',"
        f"{float(center.ra.degree):.10f},{float(center.dec.degree):.10f},{float(radius):.10f})) "
        "ORDER BY phot_g_mean_mag ASC"
    )
    params = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    url = f"https://gea.esac.esa.int/tap-server/tap/sync?{params}"
    try:
        with urllib.request.urlopen(url, timeout=int(timeout_seconds)) as response:
            payload = response.read().decode("utf-8", errors="replace")
    except Exception:
        return []
    rows = []
    for row in csv.DictReader(io.StringIO(payload)):
        try:
            rows.append(
                {
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                    "phot_g_mean_mag": float(row["phot_g_mean_mag"]),
                }
            )
        except (KeyError, TypeError, ValueError):
            continue
    return rows


def _inside_image(image: np.ndarray, x: float | None, y: float | None) -> bool:
    return bool(x is not None and y is not None and 0 <= x < image.shape[1] and 0 <= y < image.shape[0])


def _row_value(row, names: tuple[str, ...], default=None):
    columns = set(row.colnames)
    for name in names:
        if name in columns:
            return row[name]
    return default


def _plain_value(value):
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return str(value) if not isinstance(value, (str, int, float, bool, type(None))) else value


def _safe_float(value) -> float | None:
    try:
        number = float(value)
        return number if np.isfinite(number) else None
    except (TypeError, ValueError):
        return None
