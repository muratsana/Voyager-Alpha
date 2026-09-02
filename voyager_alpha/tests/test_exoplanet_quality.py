import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits

from voyager_alpha.core.exoplanet import aperture_measurement
from voyager_alpha.core.exoplanet_quality import (
    resolve_detector_linearity,
    run_photometry_preflight,
)


class ExoplanetQualityTests(unittest.TestCase):
    def test_linearity_uses_explicit_camera_limit_before_numeric_range(self):
        header = fits.Header()
        header["BITPIX"] = 16
        header["BZERO"] = 32768

        inferred = resolve_detector_linearity(header)
        explicit = resolve_detector_linearity(header, 42000)

        self.assertAlmostEqual(inferred.limit_adu, 65535 * 0.9, delta=1.0)
        self.assertFalse(inferred.verified)
        self.assertEqual(explicit.limit_adu, 42000)
        self.assertTrue(explicit.verified)

    def test_aperture_flags_a_core_above_linearity_limit(self):
        yy, xx = np.indices((80, 80), dtype=np.float32)
        image = np.full((80, 80), 1000.0, dtype=np.float32)
        image += 50000.0 * np.exp(-((xx - 40.0) ** 2 + (yy - 39.0) ** 2) / (2.0 * 1.8**2))

        measurement = aperture_measurement(
            image,
            40,
            39,
            aperture_radius=6,
            linearity_limit_adu=45000,
        )

        self.assertIn("saturated", measurement.flags)
        self.assertGreater(measurement.saturated_pixels, 0)
        self.assertGreater(measurement.fwhm_px, 2.0)

    def test_preflight_selects_stable_unsaturated_ensemble(self):
        with tempfile.TemporaryDirectory() as folder:
            files = self._write_sequence(folder, saturated_target=False)

            result = run_photometry_preflight(
                files,
                (70.0, 80.0),
                aperture_radius=6.0,
                catalog_lookup=False,
            )

            self.assertEqual(result.target.status, "PASS")
            self.assertGreaterEqual(len(result.recommended_xy), 3)
            self.assertLessEqual(len(result.recommended_xy), 8)
            self.assertTrue(all(np.hypot(x - 185.0, y - 125.0) > 3.0 for x, y in result.recommended_xy))
            self.assertTrue(any("unstable_reference" in item.flags for item in result.references))

    def test_preflight_rejects_saturated_target(self):
        with tempfile.TemporaryDirectory() as folder:
            files = self._write_sequence(folder, saturated_target=True)

            result = run_photometry_preflight(
                files,
                (70.0, 80.0),
                aperture_radius=6.0,
                catalog_lookup=False,
            )

            self.assertEqual(result.target.status, "FAIL")
            self.assertIn("saturated", result.target.flags)
            self.assertFalse(result.analysis_allowed)

    def test_preflight_recommends_spreading_an_undersampled_psf(self):
        with tempfile.TemporaryDirectory() as folder:
            files = self._write_sequence(folder, saturated_target=False, sigma=0.72)

            result = run_photometry_preflight(
                files,
                (70.0, 80.0),
                aperture_radius=5.0,
                catalog_lookup=False,
            )

            self.assertEqual(result.target.status, "CAUTION")
            self.assertIn("undersampled_psf", result.target.flags)
            self.assertTrue(any("defocus" in message for message in result.target.messages))

    @staticmethod
    def _write_sequence(folder, *, saturated_target, sigma=1.8):
        yy, xx = np.indices((190, 240), dtype=np.float32)
        stable_stars = (
            (30, 35, 10500.0),
            (118, 30, 9000.0),
            (205, 42, 8500.0),
            (35, 145, 9200.0),
            (125, 155, 10000.0),
            (210, 158, 8200.0),
        )
        files = []
        for index in range(12):
            rng = np.random.default_rng(1000 + index)
            image = rng.normal(1000.0, 5.0, yy.shape).astype(np.float32)
            target_amplitude = 35000.0 if saturated_target else 12000.0
            image += target_amplitude * np.exp(-((xx - 70.0) ** 2 + (yy - 80.0) ** 2) / (2.0 * sigma**2))
            for x, y, amplitude in stable_stars:
                image += amplitude * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * sigma**2))
            unstable_amplitude = 7000.0 * (1.0 + (0.35 if index % 2 else -0.35))
            image += unstable_amplitude * np.exp(
                -((xx - 185.0) ** 2 + (yy - 125.0) ** 2) / (2.0 * sigma**2)
            )
            path = Path(folder) / f"quality_{index:02d}.fits"
            hdu = fits.PrimaryHDU(image)
            hdu.header["DATE-OBS"] = f"2026-09-02T20:{index:02d}:00"
            hdu.header["EXPTIME"] = 60.0
            hdu.header["LINLEVEL"] = 30000.0
            hdu.header["EGAIN"] = 0.8
            hdu.writeto(path)
            files.append(str(path))
        return files


if __name__ == "__main__":
    unittest.main()
