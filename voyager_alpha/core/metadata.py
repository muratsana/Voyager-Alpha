from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from .fits_io import read_fits_header
from .models import CameraMetadata, FrameRecord


EXPOSURE_KEYS = ("EXPTIME", "EXPOSURE", "EXP_TIME")
CAMERA_KEYS = ("INSTRUME", "CAMERA", "CAMMODEL", "CCDNAME")
DETECTOR_KEYS = ("DETECTOR", "SENSOR", "SENSORID", "CCDNAME")
PIXEL_X_KEYS = ("XPIXSZ", "XPIXELSZ", "PIXSIZE1", "PIXSIZEX")
PIXEL_Y_KEYS = ("YPIXSZ", "YPIXELSZ", "PIXSIZE2", "PIXSIZEY")
BIN_X_KEYS = ("XBINNING", "CCDXBIN", "BINX", "XBIN")
BIN_Y_KEYS = ("YBINNING", "CCDYBIN", "BINY", "YBIN")
GAIN_KEYS = ("EGAIN", "GAIN", "CCDGAIN")
OFFSET_KEYS = ("OFFSET", "CCDOFFST", "BLACKLEV", "PEDESTAL")
TEMPERATURE_KEYS = ("CCD-TEMP", "CCD_TEMP", "SENSORT", "SENSOR_T", "TEMPERAT")
SET_TEMPERATURE_KEYS = ("SET-TEMP", "SET_TEMP", "CCDSET")
FOCAL_LENGTH_KEYS = ("FOCALLEN", "FOCALLENGTH", "FLENGTH")
APERTURE_KEYS = ("APTDIA", "APERTURE", "APERTURE_DIAMETER")
SCALE_KEYS = ("PIXSCALE", "SECPIX", "SCALE")


def _read_exposure(header) -> float | None:
    for key in EXPOSURE_KEYS:
        value = header.get(key)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None
    return None


def _midpoint_jd(date_obs: str | None, exposure_seconds: float | None) -> float | None:
    if not date_obs:
        return None
    try:
        start = Time(date_obs, format="fits")
        offset_days = (exposure_seconds or 0.0) / 2.0 / 86400.0
        return float(start.jd + offset_days)
    except Exception:
        return None


def _first_text(header, keys) -> str | None:
    for key in keys:
        value = header.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _first_float(header, keys) -> tuple[float | None, str | None]:
    for key in keys:
        value = header.get(key)
        if value is None:
            continue
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if number == number and abs(number) < 1e9:
            return number, key
    return None, None


def _first_int(header, keys) -> int | None:
    value, _key = _first_float(header, keys)
    if value is None or value <= 0:
        return None
    return max(1, int(round(value)))


def _read_binning(header) -> tuple[int | None, int | None]:
    bin_x = _first_int(header, BIN_X_KEYS)
    bin_y = _first_int(header, BIN_Y_KEYS)
    if bin_x is not None or bin_y is not None:
        return bin_x or bin_y, bin_y or bin_x
    text = _first_text(header, ("BINNING",))
    if text:
        normalized = text.lower().replace(" ", "").replace("*", "x")
        parts = normalized.split("x", 1)
        try:
            if len(parts) == 2:
                return max(1, int(parts[0])), max(1, int(parts[1]))
            value = max(1, int(parts[0]))
            return value, value
        except (TypeError, ValueError):
            pass
    return None, None


def _image_scale(header, pixel_size_um: float | None, binning: int | None, focal_length_mm: float | None):
    declared, key = _first_float(header, SCALE_KEYS)
    if declared is not None and 0.01 <= abs(declared) <= 100.0:
        return abs(declared), key
    try:
        wcs = WCS(header)
        if wcs.has_celestial:
            scale = float(abs(proj_plane_pixel_scales(wcs.celestial).mean()) * 3600.0)
            if 0.01 <= scale <= 100.0:
                return scale, "WCS"
    except Exception:
        pass
    if pixel_size_um and focal_length_mm and pixel_size_um > 0 and focal_length_mm > 0:
        scale = 206.265 * pixel_size_um * max(1, binning or 1) / focal_length_mm
        if 0.01 <= scale <= 100.0:
            return float(scale), "pixel/focal"
    return None, None


def read_camera_metadata(header) -> CameraMetadata:
    pixel_x, _pixel_x_key = _first_float(header, PIXEL_X_KEYS)
    pixel_y, _pixel_y_key = _first_float(header, PIXEL_Y_KEYS)
    if pixel_x is None:
        pixel_x = pixel_y
    if pixel_y is None:
        pixel_y = pixel_x
    bin_x, bin_y = _read_binning(header)
    gain, gain_keyword = _first_float(header, GAIN_KEYS)
    offset, _offset_keyword = _first_float(header, OFFSET_KEYS)
    temperature, _temperature_keyword = _first_float(header, TEMPERATURE_KEYS)
    set_temperature, _set_temperature_keyword = _first_float(header, SET_TEMPERATURE_KEYS)
    focal_length, _focal_keyword = _first_float(header, FOCAL_LENGTH_KEYS)
    aperture, _aperture_keyword = _first_float(header, APERTURE_KEYS)
    image_scale, scale_source = _image_scale(header, pixel_x, bin_x, focal_length)
    return CameraMetadata(
        instrument=_first_text(header, CAMERA_KEYS),
        detector=_first_text(header, DETECTOR_KEYS),
        pixel_size_x_um=pixel_x,
        pixel_size_y_um=pixel_y,
        binning_x=bin_x,
        binning_y=bin_y,
        gain=gain,
        gain_keyword=gain_keyword,
        offset=offset,
        sensor_temperature_c=temperature,
        set_temperature_c=set_temperature,
        filter_name=_first_text(header, ("FILTER", "FILTERID", "FILTNAME")),
        focal_length_mm=focal_length,
        aperture_mm=aperture,
        image_scale_arcsec_px=image_scale,
        image_scale_source=scale_source,
        readout_mode=_first_text(header, ("READOUTM", "READMODE", "READOUT")),
        bayer_pattern=_first_text(header, ("BAYERPAT", "BAYERPATTERN")),
    )


def inspect_frame(index: int, file_path: str) -> FrameRecord:
    header = read_fits_header(file_path)
    if int(header.get("NAXIS", 0)) != 2:
        raise ValueError(f"Only 2D FITS frames are supported: {file_path}")
    shape = (int(header["NAXIS2"]), int(header["NAXIS1"]))

    date_obs = header.get("DATE-OBS")
    exposure_seconds = _read_exposure(header)
    midpoint_jd = _midpoint_jd(date_obs, exposure_seconds)

    quality_flags = []
    if not date_obs:
        quality_flags.append("missing_DATE-OBS")
    if midpoint_jd is None:
        quality_flags.append("invalid_time")
    if exposure_seconds is None:
        quality_flags.append("missing_exposure")

    try:
        has_wcs = bool(WCS(header).has_celestial)
    except Exception:
        has_wcs = False

    if not has_wcs:
        quality_flags.append("missing_or_invalid_WCS")

    return FrameRecord(
        index=index,
        file_path=file_path,
        date_obs=date_obs,
        midpoint_jd=midpoint_jd,
        exposure_seconds=exposure_seconds,
        shape=shape,
        has_wcs=has_wcs,
        quality_flags=quality_flags,
        camera=read_camera_metadata(header),
    )


def inspect_sequence(fits_files: list[str]) -> list[FrameRecord]:
    frames = [inspect_frame(idx, file_path) for idx, file_path in enumerate(fits_files)]
    if not frames:
        return frames

    if all(frame.midpoint_jd is not None for frame in frames):
        frames = sorted(frames, key=lambda frame: frame.midpoint_jd or 0.0)
        for index, frame in enumerate(frames):
            frame.index = index

    first_shape = frames[0].shape
    for frame in frames:
        if frame.shape != first_shape:
            frame.quality_flags.append("shape_mismatch")

    exposures = [f.exposure_seconds for f in frames if f.exposure_seconds is not None]
    if exposures and max(exposures) - min(exposures) > 0.01:
        for frame in frames:
            frame.quality_flags.append("exposure_mismatch")

    return frames
