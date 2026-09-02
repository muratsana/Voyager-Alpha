from __future__ import annotations

import os

import numpy as np
from PyQt6.QtCore import QThread, pyqtSignal

from .calibration import calibrate_science_frame, load_master_dark, load_master_frame
from .exoplanet_catalog import ExoplanetCatalog
from .exoplanet import aperture_measurement, differential_light_curve_from_fluxes
from .exoplanet_quality import merge_detector_quality, run_photometry_preflight, sensor_coordinates
from .models import RegistrationSolution
from .fits_io import read_fits_image
from .metadata import inspect_sequence
from .plate_solver import AstapPlateSolver
from .registration import register_frame
from .wcs_cache import SequenceWcsResult, WcsSolutionCache, propagate_wcs_header


class CatalogUpdateWorker(QThread):
    progress = pyqtSignal(int, str, str)
    log = pyqtSignal(str)
    result_ready = pyqtSignal(object)
    failed = pyqtSignal(str)
    finished_scan = pyqtSignal()

    def __init__(self, catalog: ExoplanetCatalog, parent=None):
        super().__init__(parent)
        self.catalog = catalog
        self.is_running = True

    def run(self):
        try:
            counts = self.catalog.update_all(
                progress_callback=self.progress.emit,
                stop_callback=lambda: not self.is_running,
            )
            self.result_ready.emit(counts)
            self.log.emit(f"INFO|Catalog|Transit catalogs updated: {sum(counts.values()):,} source records")
        except Exception as exc:
            self.log.emit(f"ERROR|Catalog|{exc}")
            self.failed.emit(str(exc))
        finally:
            self.finished_scan.emit()

    def stop(self):
        self.is_running = False


class SequencePlateSolveWorker(QThread):
    progress = pyqtSignal(int, str)
    log = pyqtSignal(str)
    result_ready = pyqtSignal(object)
    finished_scan = pyqtSignal()

    def __init__(self, fits_files, *, cache=None, solver=None, parent=None):
        super().__init__(parent)
        self.fits_files = list(fits_files)
        self.cache = cache or WcsSolutionCache()
        self.solver = solver or AstapPlateSolver()
        self.is_running = True

    def run(self):
        result = SequenceWcsResult()
        try:
            if not self.fits_files:
                raise ValueError("Plate solve için FITS karesi bulunamadı.")
            for file_path in self.fits_files:
                header, method = self.cache.resolve(file_path)
                if header is not None:
                    result.headers[file_path] = header
                    result.methods[file_path] = method
            if len(result.headers) == len(self.fits_files):
                self.progress.emit(100, "Tüm karelerin WCS çözümü hazır")
                self.result_ready.emit(result)
                return
            if not self.solver.is_available():
                raise RuntimeError("ASTAP bulunamadı; WCS'siz sekans çözülemedi.")

            reference_path = next((path for path in self.fits_files if path in result.headers), self.fits_files[0])
            if reference_path not in result.headers:
                self.progress.emit(2, "Referans kare ASTAP ile çözülüyor")
                solved = self.solver.solve(reference_path)
                if not solved.success or solved.header is None:
                    raise RuntimeError(solved.message or "Referans kare ASTAP ile çözülemedi.")
                self.cache.store(reference_path, solved.header)
                result.headers[reference_path] = solved.header
                result.methods[reference_path] = "astap"
                self.log.emit(f"INFO|Astrometry|Reference solved with ASTAP: {os.path.basename(reference_path)}")

            reference_image = read_fits_image(reference_path)
            reference_header = result.headers[reference_path]
            total = len(self.fits_files)
            for index, file_path in enumerate(self.fits_files, start=1):
                if not self.is_running:
                    return
                if file_path in result.headers:
                    self.progress.emit(int(index / total * 100), f"WCS hazır {index}/{total}")
                    continue
                try:
                    image = read_fits_image(file_path)
                    _aligned, registration = register_frame(reference_image, image)
                    reliable = registration.matched_stars >= 8 and registration.rms_px <= 2.0
                    if reliable:
                        header = propagate_wcs_header(reference_header, registration)
                        method = f"propagated:{registration.rms_px:.2f}px"
                    else:
                        solved = self.solver.solve(file_path)
                        if not solved.success or solved.header is None:
                            raise RuntimeError(solved.message or "ASTAP çözüm üretemedi.")
                        header = solved.header
                        method = "astap"
                    self.cache.store(file_path, header)
                    result.headers[file_path] = header
                    result.methods[file_path] = method
                    level = "INFO" if method.startswith("propagated") else "WARN"
                    self.log.emit(f"{level}|Astrometry|{os.path.basename(file_path)}: {method}")
                except Exception as exc:
                    result.failures[file_path] = str(exc)
                    self.log.emit(f"ERROR|Astrometry|{os.path.basename(file_path)}: {exc}")
                self.progress.emit(int(index / total * 100), f"Astrometry {index}/{total}")
            self.result_ready.emit(result)
        except Exception as exc:
            for file_path in self.fits_files:
                if file_path not in result.headers:
                    result.failures.setdefault(file_path, str(exc))
            self.log.emit(f"ERROR|Astrometry|{exc}")
            self.result_ready.emit(result)
        finally:
            self.finished_scan.emit()

    def stop(self):
        self.is_running = False


class PhotometryPreflightWorker(QThread):
    progress = pyqtSignal(int, str)
    log = pyqtSignal(str)
    result_ready = pyqtSignal(object)
    failed = pyqtSignal(str)
    finished_scan = pyqtSignal()

    def __init__(
        self,
        fits_files,
        target_xy,
        *,
        aperture_radius=6.0,
        linearity_limit_adu=0.0,
        master_paths=None,
        reference_header=None,
        catalog_lookup=True,
        parent=None,
    ):
        super().__init__(parent)
        self.fits_files = list(fits_files)
        self.target_xy = tuple(target_xy)
        self.aperture_radius = float(aperture_radius)
        self.linearity_limit_adu = float(linearity_limit_adu)
        self.master_paths = dict(master_paths or {})
        self.reference_header = reference_header
        self.catalog_lookup = bool(catalog_lookup)
        self.is_running = True

    def run(self):
        try:
            result = run_photometry_preflight(
                self.fits_files,
                self.target_xy,
                aperture_radius=self.aperture_radius,
                linearity_limit_adu=self.linearity_limit_adu,
                master_paths=self.master_paths,
                reference_header=self.reference_header,
                catalog_lookup=self.catalog_lookup,
                progress_callback=self.progress.emit,
                stop_callback=lambda: not self.is_running,
            )
            if self.is_running:
                self.result_ready.emit(result)
                self.log.emit(
                    f"INFO|Photometry QC|Target {result.target.status}; "
                    f"{len(result.recommended_xy)} comparison stars selected"
                )
        except Exception as exc:
            if self.is_running:
                self.failed.emit(str(exc))
                self.log.emit(f"ERROR|Photometry QC|{exc}")
        finally:
            self.finished_scan.emit()

    def stop(self):
        self.is_running = False


class ExoplanetWorker(QThread):
    progress = pyqtSignal(int, str)
    file_progress = pyqtSignal(int, int, str, str)
    log = pyqtSignal(str)
    measurement_ready = pyqtSignal(int, float, float)
    result_ready = pyqtSignal(object)
    finished_scan = pyqtSignal()

    def __init__(
        self,
        fits_files,
        target_xy,
        comparison_xy,
        *,
        aperture_radius=6.0,
        detrend_order=2,
        recenter=True,
        linearity_limit_adu=None,
        master_paths=None,
        parent=None,
    ):
        super().__init__(parent)
        self.fits_files = list(fits_files)
        self.target_xy = tuple(target_xy)
        self.comparison_xy = [tuple(value) for value in comparison_xy]
        self.aperture_radius = float(aperture_radius)
        self.detrend_order = int(detrend_order)
        self.recenter = bool(recenter)
        self.linearity_limit_adu = (
            float(linearity_limit_adu) if linearity_limit_adu is not None and linearity_limit_adu > 0 else None
        )
        self.master_paths = dict(master_paths or {})
        self.is_running = True

    def run(self):
        try:
            records = inspect_sequence(self.fits_files)
            if len(records) < 5:
                raise ValueError("Transit analizi için en az 5 FITS karesi gerekir.")
            if any(record.midpoint_jd is None for record in records):
                raise ValueError("Tüm karelerde geçerli DATE-OBS ve poz süresi bulunmalıdır.")
            if not self.comparison_xy:
                raise ValueError("En az bir karşılaştırma yıldızı seçilmelidir.")
            records.sort(key=lambda record: record.midpoint_jd or 0.0)
            files = [record.file_path for record in records]
            total = len(files)

            bias = load_master_frame(self.master_paths["bias"]) if self.master_paths.get("bias") else None
            flat = load_master_frame(self.master_paths["flat"]) if self.master_paths.get("flat") else None
            dark = None
            dark_exposure = None
            if self.master_paths.get("dark"):
                dark, dark_exposure = load_master_dark(self.master_paths["dark"])
            selected_masters = [name for name in ("bias", "dark", "flat") if self.master_paths.get(name)]
            if selected_masters:
                self.log.emit(f"INFO|Calibration|Applied masters: {', '.join(selected_masters)}")

            target_values = []
            target_uncertainties = []
            comparison_values = []
            comparison_uncertainties = []
            measurement_flags = set()
            centroid_offsets = []
            saturated_target_frames = []
            reference = None
            for index, (file_path, record) in enumerate(zip(files, records)):
                if not self.is_running:
                    return
                self.file_progress.emit(index + 1, total, os.path.basename(file_path), "Calibrate and register")
                raw = np.asarray(read_fits_image(file_path), dtype=np.float32)
                data = calibrate_science_frame(
                    raw,
                    master_bias=bias,
                    master_dark=dark,
                    master_flat=flat,
                    science_exposure=record.exposure_seconds,
                    dark_exposure=dark_exposure,
                )
                if reference is None:
                    reference = data
                    aligned = reference
                    solution = RegistrationSolution.identity()
                    self.log.emit("INFO|Align|Reference frame established")
                else:
                    aligned, solution = register_frame(reference, data)
                    level = "INFO" if solution.rms_px <= 1.5 else "WARN"
                    self.log.emit(
                        f"{level}|Align|Frame {index + 1}: {solution.method}, "
                        f"RMS {solution.rms_px:.2f} px, {solution.matched_stars} stars"
                    )
                target = aperture_measurement(
                    aligned,
                    *self.target_xy,
                    aperture_radius=self.aperture_radius,
                    recenter=self.recenter,
                    gain_e_per_adu=_usable_gain(record),
                )
                target_sensor_xy = sensor_coordinates(self.target_xy, solution)
                target = merge_detector_quality(
                    target,
                    aperture_measurement(
                        raw,
                        *target_sensor_xy,
                        aperture_radius=self.aperture_radius,
                        recenter=self.recenter,
                        linearity_limit_adu=self.linearity_limit_adu,
                    ),
                )
                comparisons = [
                    merge_detector_quality(
                        aperture_measurement(
                            aligned,
                            *xy,
                            aperture_radius=self.aperture_radius,
                            recenter=self.recenter,
                            gain_e_per_adu=_usable_gain(record),
                        ),
                        aperture_measurement(
                            raw,
                            *sensor_coordinates(xy, solution),
                            aperture_radius=self.aperture_radius,
                            recenter=self.recenter,
                            linearity_limit_adu=self.linearity_limit_adu,
                        ),
                    )
                    for xy in self.comparison_xy
                ]
                if "saturated" in target.flags:
                    saturated_target_frames.append(index + 1)
                target_values.append(target.flux)
                target_uncertainties.append(target.uncertainty)
                comparison_values.append(
                    [float("nan") if "saturated" in value.flags else value.flux for value in comparisons]
                )
                comparison_uncertainties.append(
                    [float("nan") if "saturated" in value.flags else value.uncertainty for value in comparisons]
                )
                measurement_flags.update(target.flags)
                for value in comparisons:
                    measurement_flags.update(value.flags)
                centroid_offsets.append(float(np.hypot(target.x - self.target_xy[0], target.y - self.target_xy[1])))
                comparison_sum = float(np.nansum([value.flux for value in comparisons]))
                self.measurement_ready.emit(index, float(target.flux), comparison_sum)
                self.progress.emit(8 + int(((index + 1) / total) * 82), "Photometry in progress")

            if not self.is_running:
                return
            if saturated_target_frames:
                raise ValueError(
                    "Hedef yıldız kamera lineerlik sınırını aşıyor: kare "
                    + ", ".join(map(str, saturated_target_frames[:12]))
                    + ("..." if len(saturated_target_frames) > 12 else "")
                )
            self.progress.emit(94, "Differential light curve is being calculated")
            result = differential_light_curve_from_fluxes(
                [record.midpoint_jd for record in records],
                target_values,
                comparison_values,
                target_uncertainties=target_uncertainties,
                comparison_uncertainties=comparison_uncertainties,
                detrend_order=self.detrend_order,
            )
            result.quality_flags.extend(sorted(measurement_flags))
            if centroid_offsets:
                self.log.emit(
                    f"INFO|Centroid|Target recenter median {np.median(centroid_offsets):.2f} px, "
                    f"max {np.max(centroid_offsets):.2f} px"
                )
            self.result_ready.emit(result)
            level = "WARN" if result.transit_candidate else "INFO"
            self.log.emit(f"{level}|Transit|{result.message}")
            self.progress.emit(100, "Transit analysis complete")
        except Exception as exc:
            self.log.emit(f"ERROR|Transit|{exc}")
        finally:
            self.finished_scan.emit()

    def stop(self):
        self.is_running = False


def _usable_gain(record):
    gain = record.camera.gain
    if gain is None or gain <= 0:
        return None
    if record.camera.gain_keyword == "EGAIN" or gain <= 10.0:
        return float(gain)
    return None
