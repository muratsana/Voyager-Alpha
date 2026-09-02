from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord, get_sun
from astropy.time import Time


SUPPORTED_TIME_SCALES = {"utc", "tai", "tt", "tdb", "tcg", "tcb", "ut1"}


@dataclass
class FrameTimeInfo:
    midpoint: Time | None
    source: str | None
    declared_scale: str
    flags: list[str] = field(default_factory=list)
    header_bjd_tdb: float | None = None

    @property
    def jd_utc(self) -> float | None:
        return float(self.midpoint.utc.jd) if self.midpoint is not None else None

    @property
    def midpoint_utc(self) -> str | None:
        if self.midpoint is None:
            return None
        value = self.midpoint.utc
        value.precision = 3
        return f"{value.isot}Z"


@dataclass
class ObservingGeometry:
    bjd_tdb: float | None = None
    airmass: float | None = None
    altitude_deg: float | None = None
    sun_altitude_deg: float | None = None
    flags: list[str] = field(default_factory=list)


def parse_frame_time(header, exposure_seconds: float | None = None) -> FrameTimeInfo:
    """Resolve a FITS exposure midpoint without silently changing time scales."""

    raw_scale = str(header.get("TIMESYS", "UTC") or "UTC").strip().lower()
    flags: list[str] = []
    if raw_scale not in SUPPORTED_TIME_SCALES:
        flags.append(f"unsupported_TIMESYS_{raw_scale}")
        scale = "utc"
    else:
        scale = raw_scale

    header_bjd = _first_numeric(header, ("BJD_TDB", "BJD-TDB"))
    direct_candidates = (
        ("JD-AVG", "jd"),
        ("JD-MID", "jd"),
        ("MJD-AVG", "mjd"),
        ("MJD-MID", "mjd"),
        ("DATE-AVG", None),
        ("DATE-MID", None),
    )
    for key, time_format in direct_candidates:
        if header.get(key) is None:
            continue
        midpoint = _parse_time(header.get(key), scale, time_format)
        if midpoint is not None:
            return FrameTimeInfo(midpoint, key, scale.upper(), flags, header_bjd)

    start_key = next((key for key in ("DATE-BEG", "DATE-OBS") if header.get(key) is not None), None)
    end_key = next((key for key in ("DATE-END",) if header.get(key) is not None), None)
    if start_key and end_key:
        start = _parse_time(header.get(start_key), scale)
        end = _parse_time(header.get(end_key), scale)
        if start is not None and end is not None and end > start:
            return FrameTimeInfo(start + (end - start) / 2.0, f"{start_key}+{end_key}", scale.upper(), flags, header_bjd)

    numeric_start = (
        ("JD-OBS", "jd"),
        ("MJD-OBS", "mjd"),
        ("JD", "jd"),
        ("MJD", "mjd"),
    )
    for key, time_format in numeric_start:
        if header.get(key) is None:
            continue
        start = _parse_time(header.get(key), scale, time_format)
        if start is not None:
            if key in {"JD", "MJD"}:
                flags.append(f"ambiguous_{key}_assumed_start")
            return FrameTimeInfo(
                _offset_to_midpoint(start, exposure_seconds),
                key,
                scale.upper(),
                flags,
                header_bjd,
            )

    if start_key:
        start = _parse_time(header.get(start_key), scale)
        if start is not None:
            if start_key == "DATE-OBS":
                flags.append("date_obs_assumed_exposure_start")
            return FrameTimeInfo(
                _offset_to_midpoint(start, exposure_seconds),
                start_key,
                scale.upper(),
                flags,
                header_bjd,
            )

    if header.get("DATE-END") is not None:
        end = _parse_time(header.get("DATE-END"), scale)
        if end is not None:
            seconds = max(float(exposure_seconds or 0.0), 0.0)
            return FrameTimeInfo(end - seconds * u.s / 2.0, "DATE-END", scale.upper(), flags, header_bjd)

    flags.append("invalid_time")
    return FrameTimeInfo(None, None, scale.upper(), flags, header_bjd)


def observatory_location_from_header(header) -> EarthLocation | None:
    latitude = _first_angle(
        header,
        ("SITELAT", "SITELATITUDE", "LAT-OBS", "OBSLAT", "OBS-LAT", "LATITUDE"),
        latitude=True,
    )
    longitude = _first_angle(
        header,
        ("SITELONG", "SITELON", "LONG-OBS", "OBSLONG", "OBS-LONG", "LONGITUD"),
        latitude=False,
    )
    if latitude is None or longitude is None:
        return None
    height = _first_numeric(header, ("SITEELEV", "SITEALT", "ALT-OBS", "OBSALT", "ELEVATIO"))
    return EarthLocation.from_geodetic(longitude, latitude, float(height or 0.0) * u.m)


def target_coord_from_header(header) -> SkyCoord | None:
    ra_value = next((header.get(key) for key in ("OBJCTRA", "OBJRA", "RA") if header.get(key) is not None), None)
    dec_value = next((header.get(key) for key in ("OBJCTDEC", "OBJDEC", "DEC") if header.get(key) is not None), None)
    if ra_value is None or dec_value is None:
        return None
    try:
        if ":" in str(ra_value):
            return SkyCoord(ra_value, dec_value, unit=(u.hourangle, u.deg), frame="icrs")
        return SkyCoord(float(ra_value) * u.deg, float(dec_value) * u.deg, frame="icrs")
    except (TypeError, ValueError):
        return None


def observing_geometry(
    midpoint: Time | None,
    target: SkyCoord | None,
    location: EarthLocation | None,
) -> ObservingGeometry:
    if midpoint is None:
        return ObservingGeometry(flags=["missing_midpoint_time"])
    if target is None:
        return ObservingGeometry(flags=["missing_target_coordinates"])
    if location is None:
        return ObservingGeometry(flags=["missing_observatory_location"])

    time_at_site = Time(midpoint.utc.jd, format="jd", scale="utc", location=location)
    light_time = time_at_site.light_travel_time(target, kind="barycentric")
    bjd_tdb = float((time_at_site.tdb + light_time).jd)
    horizontal = target.transform_to(AltAz(obstime=time_at_site, location=location, pressure=0 * u.hPa))
    altitude = float(horizontal.alt.deg)
    secz = float(horizontal.secz.value)
    airmass = secz if altitude > 0.0 and np.isfinite(secz) and secz >= 1.0 else None
    sun_altitude = float(get_sun(time_at_site).transform_to(AltAz(obstime=time_at_site, location=location)).alt.deg)
    flags = [] if airmass is not None else ["target_below_horizon"]
    return ObservingGeometry(bjd_tdb, airmass, altitude, sun_altitude, flags)


def _offset_to_midpoint(start: Time, exposure_seconds: float | None) -> Time:
    seconds = max(float(exposure_seconds or 0.0), 0.0)
    return start + seconds * u.s / 2.0


def _parse_time(value, scale: str, time_format: str | None = None) -> Time | None:
    try:
        if time_format is not None:
            return Time(float(value), format=time_format, scale=scale)
        text = str(value).strip()
        try:
            return Time(text, format="fits", scale=scale)
        except (TypeError, ValueError):
            return Time(text, scale=scale)
    except (TypeError, ValueError, OverflowError):
        return None


def _first_numeric(header, keys) -> float | None:
    for key in keys:
        value = header.get(key)
        if value is None:
            continue
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(number):
            return number
    return None


def _first_angle(header, keys, *, latitude: bool) -> Angle | None:
    for key in keys:
        value = header.get(key)
        if value is None:
            continue
        try:
            if isinstance(value, str) and ":" in value:
                angle = Angle(value, unit=u.deg)
            else:
                angle = Angle(float(value), unit=u.deg)
        except (TypeError, ValueError, u.UnitsError):
            continue
        degrees = float(angle.deg)
        if latitude and not -90.0 <= degrees <= 90.0:
            continue
        if not latitude and not -360.0 <= degrees <= 360.0:
            continue
        return angle
    return None
