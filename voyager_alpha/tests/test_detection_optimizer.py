import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from voyager_alpha.core.detection_optimizer import optimize_detection_settings
from voyager_alpha.core.pipeline import AsteroidWorker


class DetectionOptimizerTests(unittest.TestCase):
    def test_builds_image_aware_profile_from_timed_wcs_sequence(self):
        with tempfile.TemporaryDirectory() as folder:
            files = []
            for index in range(7):
                image = self._frame(index)
                path = Path(folder) / f"frame_{index:03d}.fits"
                fits.writeto(path, image, self._header(image.shape, index), overwrite=True)
                files.append(str(path))

            result = optimize_detection_settings(files)

            self.assertEqual(result.metrics["frame_count"], 7)
            self.assertAlmostEqual(float(result.metrics["pixel_scale_arcsec"]), 1.2, places=1)
            self.assertGreater(float(result.metrics["fwhm_px"]), 2.0)
            self.assertGreaterEqual(float(result.settings["sigma"]), 4.5)
            self.assertLessEqual(float(result.settings["sigma"]), 7.0)
            self.assertEqual(result.settings["min_frames"], 4)
            self.assertLess(float(result.settings["max_fit_rms_px"]), 2.0)
            self.assertGreater(float(result.settings["max_step_px"]), float(result.settings["min_motion_px_per_frame"]))

            settings = dict(result.settings)
            pipeline_results = []
            worker = AsteroidWorker(
                files,
                sigma=settings.pop("sigma"),
                min_pix=settings.pop("min_pixels"),
                min_tracklet_frames=settings.pop("min_frames"),
                max_sources_per_frame=settings.pop("max_sources"),
                match_known_objects=False,
                **settings,
            )
            worker.result_ready.connect(pipeline_results.append)
            worker.run()
            self.assertEqual(len(pipeline_results), 1)
            self.assertEqual(len(pipeline_results[0].tracklets), 1)
            self.assertEqual(pipeline_results[0].tracklets[0].frames_detected, 7)

    def _frame(self, index):
        rng = np.random.default_rng(900 + index)
        image = rng.normal(1000.0, 4.0, (160, 180)).astype(np.float32)
        for x, y, amplitude in ((30, 35, 500), (70, 120, 700), (120, 55, 620), (150, 135, 540)):
            self._add_gaussian(image, x, y, amplitude)
        self._add_gaussian(image, 45 + 2.2 * index, 65 + 1.1 * index, 90)
        return image

    @staticmethod
    def _add_gaussian(image, x, y, amplitude):
        yy, xx = np.indices(image.shape)
        image += amplitude * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.35**2))

    @staticmethod
    def _header(shape, index):
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [shape[1] / 2, shape[0] / 2]
        wcs.wcs.cdelt = np.asarray([-1.2 / 3600.0, 1.2 / 3600.0])
        wcs.wcs.crval = [150.0, 22.0]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        header = wcs.to_header()
        header["DATE-OBS"] = f"2026-09-01T20:0{index}:00"
        header["EXPTIME"] = 60.0
        return header


if __name__ == "__main__":
    unittest.main()
