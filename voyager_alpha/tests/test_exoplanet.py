import unittest
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits

from voyager_alpha.core.exoplanet import (
    aperture_measurement,
    differential_light_curve,
    differential_light_curve_from_fluxes,
    fit_limb_darkened_transit,
)
from voyager_alpha.core.exoplanet_worker import ExoplanetWorker


class ExoplanetTests(unittest.TestCase):
    def test_detects_synthetic_transit(self):
        frames = []
        times = []
        transit_scale = [1.0] * 5 + [0.94] * 5 + [1.0] * 5
        for index, scale in enumerate(transit_scale):
            image = np.full((80, 80), 100.0, dtype=np.float32)
            self._add_star(image, 30, 40, 6000.0 * scale)
            self._add_star(image, 55, 35, 5000.0)
            frames.append(image)
            times.append(2460000.0 + index / 1440.0)

        result = differential_light_curve(frames, times, (30, 40), [(55, 35)], aperture_radius=5.0)

        self.assertTrue(result.transit_candidate)
        self.assertGreater(result.depth, 0.01)

    def test_streaming_worker_returns_light_curve(self):
        with tempfile.TemporaryDirectory() as folder:
            files = []
            for index, scale in enumerate([1.0] * 5 + [0.94] * 5 + [1.0] * 5):
                image = np.full((80, 80), 100.0, dtype=np.float32)
                self._add_star(image, 30, 40, 6000.0 * scale)
                self._add_star(image, 55, 35, 5000.0)
                path = Path(folder) / f"frame_{index:02d}.fits"
                hdu = fits.PrimaryHDU(image)
                hdu.header["DATE-OBS"] = f"2026-09-01T20:{index:02d}:00"
                hdu.header["EXPTIME"] = 60.0
                hdu.writeto(path)
                files.append(str(path))

            results = []
            worker = ExoplanetWorker(files, (30, 40), [(55, 35)], aperture_radius=5.0)
            worker.result_ready.connect(results.append)
            worker.run()

            self.assertEqual(len(results), 1)
            self.assertTrue(results[0].transit_candidate)

    def test_streaming_worker_blocks_a_saturated_target(self):
        with tempfile.TemporaryDirectory() as folder:
            files = []
            for index in range(8):
                image = np.full((80, 80), 100.0, dtype=np.float32)
                self._add_star(image, 30, 40, 50000.0)
                self._add_star(image, 55, 35, 9000.0)
                path = Path(folder) / f"saturated_{index:02d}.fits"
                hdu = fits.PrimaryHDU(image)
                hdu.header["DATE-OBS"] = f"2026-09-01T21:{index:02d}:00"
                hdu.header["EXPTIME"] = 60.0
                hdu.writeto(path)
                files.append(str(path))

            results = []
            logs = []
            worker = ExoplanetWorker(
                files,
                (30, 40),
                [(55, 35)],
                aperture_radius=5.0,
                linearity_limit_adu=45000.0,
            )
            worker.result_ready.connect(results.append)
            worker.log.connect(logs.append)
            worker.run()

            self.assertEqual(results, [])
            self.assertTrue(any("lineerlik sınırını aşıyor" in message for message in logs))

    def test_unstable_comparison_star_is_downweighted(self):
        times = [2460000.0 + index / 1440.0 for index in range(15)]
        target = np.asarray([10000.0] * 5 + [9600.0] * 5 + [10000.0] * 5)
        stable = np.full(15, 8000.0)
        unstable = np.asarray([7000.0, 9000.0] * 7 + [7000.0])

        result = differential_light_curve_from_fluxes(
            times,
            target,
            np.column_stack((stable, unstable)),
            target_uncertainties=np.full(15, 20.0),
            comparison_uncertainties=np.full((15, 2), 20.0),
        )

        self.assertTrue(result.transit_candidate)
        self.assertGreater(result.comparison_weights[0], result.comparison_weights[1] * 10)
        self.assertIn("comparison_2_unstable", result.quality_flags)

    def test_limb_darkened_fit_recovers_midpoint(self):
        times = np.linspace(2460000.0, 2460000.03, 120)
        midpoint = 2460000.015
        minutes = (times - midpoint) * 1440.0
        contact = np.sqrt((1.0 + 0.1) ** 2 - 0.25**2)
        velocity = 2.0 * contact / 18.0
        separation = np.sqrt(0.25**2 + (minutes * velocity) ** 2)
        from exoplanet_core import quad_limbdark_light_curve

        flux = 1.0 + quad_limbdark_light_curve(0.3, 0.2, separation, 0.1)
        uncertainty = np.full(len(times), 0.0005)
        fit = fit_limb_darkened_transit(
            times,
            flux,
            uncertainty,
            initial_midpoint=midpoint,
            initial_duration_minutes=18.0,
            initial_depth=0.01,
        )

        self.assertTrue(fit.success)
        self.assertAlmostEqual(fit.mid_transit_jd, midpoint, delta=1.0 / 1440.0)
        self.assertAlmostEqual(fit.radius_ratio, 0.1, delta=0.03)

    def test_flat_light_curve_is_not_a_candidate(self):
        rng = np.random.default_rng(123)
        times = [2460000.0 + index / 1440.0 for index in range(30)]
        target = 10000.0 + rng.normal(0.0, 8.0, 30)
        comparison = 8000.0 + rng.normal(0.0, 7.0, 30)

        result = differential_light_curve_from_fluxes(
            times,
            target,
            comparison,
            target_uncertainties=np.full(30, 10.0),
            comparison_uncertainties=np.full(30, 10.0),
        )

        self.assertFalse(result.transit_candidate)

    def test_aperture_measurement_recenters_star(self):
        yy, xx = np.indices((70, 70), dtype=np.float32)
        image = np.full((70, 70), 100.0, dtype=np.float32)
        image += 5000.0 * np.exp(-((xx - 31.2) ** 2 + (yy - 39.4) ** 2) / (2 * 1.7**2))

        measurement = aperture_measurement(image, 30.0, 40.0, aperture_radius=5.0, recenter=True)

        self.assertAlmostEqual(measurement.x, 31.2, delta=0.35)
        self.assertAlmostEqual(measurement.y, 39.4, delta=0.35)
        self.assertGreater(measurement.flux, 0.0)

    def _add_star(self, image, x, y, amplitude):
        yy, xx = np.indices(image.shape)
        image += amplitude * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 1.8**2))
