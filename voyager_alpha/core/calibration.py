from __future__ import annotations

import numpy as np

from .fits_io import read_fits_image


def load_master_dark(file_path: str) -> tuple[np.ndarray, float | None]:
    data, header = read_fits_image(file_path, header=True)
    exposure = header.get("EXPTIME", header.get("EXPOSURE"))
    try:
        exposure_seconds = float(exposure) if exposure is not None else None
    except (TypeError, ValueError):
        exposure_seconds = None
    return np.asarray(data, dtype=np.float32), exposure_seconds


def load_master_frame(file_path: str) -> np.ndarray:
    return np.asarray(read_fits_image(file_path), dtype=np.float32)


def calibrate_with_master_dark(
    science: np.ndarray,
    master_dark: np.ndarray | None,
    *,
    science_exposure: float | None = None,
    dark_exposure: float | None = None,
) -> np.ndarray:
    calibrated = np.asarray(science, dtype=np.float32)
    if master_dark is None:
        return calibrated.copy()
    dark = np.asarray(master_dark, dtype=np.float32)
    if dark.shape != calibrated.shape:
        raise ValueError(
            f"Master dark boyutu science frame ile uyuşmuyor: {dark.shape} != {calibrated.shape}"
        )
    scale = 1.0
    if science_exposure and dark_exposure and dark_exposure > 0:
        scale = float(science_exposure) / float(dark_exposure)
    return calibrated - dark * scale


def calibrate_science_frame(
    science: np.ndarray,
    *,
    master_bias: np.ndarray | None = None,
    master_dark: np.ndarray | None = None,
    master_flat: np.ndarray | None = None,
    science_exposure: float | None = None,
    dark_exposure: float | None = None,
) -> np.ndarray:
    """Apply bias, exposure-scaled dark and normalized flat calibration."""

    calibrated = np.asarray(science, dtype=np.float32).copy()
    for name, master in (
        ("bias", master_bias),
        ("dark", master_dark),
        ("flat", master_flat),
    ):
        if master is not None and np.asarray(master).shape != calibrated.shape:
            raise ValueError(
                f"Master {name} boyutu science frame ile uyuşmuyor: "
                f"{np.asarray(master).shape} != {calibrated.shape}"
            )

    if master_bias is not None:
        calibrated -= np.asarray(master_bias, dtype=np.float32)

    if master_dark is not None:
        dark = np.asarray(master_dark, dtype=np.float32)
        scale = 1.0
        if science_exposure and dark_exposure and dark_exposure > 0:
            scale = float(science_exposure) / float(dark_exposure)
        calibrated -= dark * scale

    if master_flat is not None:
        flat = np.asarray(master_flat, dtype=np.float32)
        finite_positive = flat[np.isfinite(flat) & (flat > 0)]
        if finite_positive.size == 0:
            raise ValueError("Master flat geçerli pozitif piksel içermiyor.")
        normalization = float(np.nanmedian(finite_positive))
        normalized_flat = flat / max(normalization, 1e-12)
        valid = np.isfinite(normalized_flat) & (normalized_flat > 0.05)
        fill = float(np.nanmedian(calibrated[np.isfinite(calibrated)]))
        calibrated = np.divide(
            calibrated,
            normalized_flat,
            out=np.full_like(calibrated, fill, dtype=np.float32),
            where=valid,
        )

    return np.nan_to_num(calibrated, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
