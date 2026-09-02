from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .fits_io import read_fits_header
from .models import RegistrationSolution


@dataclass
class SequenceWcsResult:
    headers: dict[str, fits.Header] = field(default_factory=dict)
    methods: dict[str, str] = field(default_factory=dict)
    failures: dict[str, str] = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        return bool(self.headers) and not self.failures


class WcsSolutionCache:
    def __init__(self, cache_directory: str | Path | None = None):
        self.cache_directory = Path(cache_directory) if cache_directory else default_wcs_cache_directory()
        self.cache_directory.mkdir(parents=True, exist_ok=True)

    def load(self, fits_path: str) -> fits.Header | None:
        path = self._path_for(fits_path)
        if not path.is_file():
            return None
        try:
            header = fits.Header.fromstring(path.read_text(encoding="ascii"), sep="\n")
            return header if has_celestial_wcs(header) else None
        except Exception:
            return None

    def store(self, fits_path: str, header: fits.Header):
        if not has_celestial_wcs(header):
            raise ValueError("Geçerli göksel WCS cache'e yazılamaz.")
        target = self._path_for(fits_path)
        temporary = target.with_suffix(".tmp")
        temporary.write_text(header.tostring(sep="\n", endcard=True, padding=False), encoding="ascii")
        temporary.replace(target)

    def resolve(self, fits_path: str) -> tuple[fits.Header | None, str]:
        try:
            header = read_fits_header(fits_path)
            if has_celestial_wcs(header):
                return header, "fits-header"
        except Exception:
            pass
        cached = self.load(fits_path)
        return (cached, "cache") if cached is not None else (None, "missing")

    def _path_for(self, fits_path: str) -> Path:
        source = Path(fits_path).resolve()
        stat = source.stat()
        signature = f"{source}|{stat.st_size}|{stat.st_mtime_ns}".encode("utf-8", errors="surrogatepass")
        digest = hashlib.sha256(signature).hexdigest()[:32]
        return self.cache_directory / f"{digest}.wcs"


def default_wcs_cache_directory() -> Path:
    local = os.environ.get("LOCALAPPDATA")
    root = Path(local) if local else Path.home() / ".voyager-alpha"
    return root / "Astrohub" / "VoyagerAlpha" / "wcs"


def has_celestial_wcs(header: fits.Header | None) -> bool:
    if header is None:
        return False
    try:
        return bool(WCS(header).has_celestial)
    except Exception:
        return False


def propagate_wcs_header(
    reference_header: fits.Header,
    registration: RegistrationSolution,
) -> fits.Header:
    reference = WCS(reference_header).celestial
    if not reference.has_celestial:
        raise ValueError("Referans başlıkta göksel WCS bulunamadı.")
    matrix = np.asarray(registration.matrix_xy, dtype=np.float64)
    offset = np.asarray(registration.offset_xy, dtype=np.float64)
    if matrix.shape != (2, 2) or offset.shape != (2,) or abs(np.linalg.det(matrix)) < 1e-8:
        raise ValueError("WCS taşımak için geçersiz kayıt dönüşümü.")

    reference_crpix_zero = np.asarray(reference.wcs.crpix, dtype=np.float64) - 1.0
    moving_crpix_zero = np.linalg.solve(matrix, reference_crpix_zero - offset)
    linear = np.asarray(reference.pixel_scale_matrix, dtype=np.float64) @ matrix

    moving = WCS(naxis=2)
    moving.wcs.crpix = moving_crpix_zero + 1.0
    moving.wcs.crval = np.asarray(reference.wcs.crval, dtype=np.float64)
    moving.wcs.ctype = list(reference.wcs.ctype)
    moving.wcs.cunit = list(reference.wcs.cunit)
    moving.wcs.cd = linear
    if reference.wcs.radesys:
        moving.wcs.radesys = reference.wcs.radesys
    if np.isfinite(reference.wcs.equinox):
        moving.wcs.equinox = reference.wcs.equinox
    return moving.to_header(relax=True)


def pixel_to_sky(header: fits.Header, x: float, y: float) -> tuple[float, float]:
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError("Göksel WCS bulunamadı.")
    ra, dec = wcs.all_pix2world(float(x), float(y), 0)
    return float(ra), float(dec)
