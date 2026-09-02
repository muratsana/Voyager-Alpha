import unittest

import numpy as np

from voyager_alpha.core.detection import extract_residual_detections, suppress_persistent_residuals
from voyager_alpha.core.models import Detection


class DetectionArtifactTests(unittest.TestCase):
    def test_rejects_single_pixel_and_compact_plateau_but_keeps_psf(self):
        rng = np.random.default_rng(24)
        pixel = rng.normal(0.0, 1.0, (80, 80)).astype(np.float32)
        pixel[40, 40] += 60.0
        block = rng.normal(0.0, 1.0, (80, 80)).astype(np.float32)
        block[40:42, 40:42] += 30.0
        source = rng.normal(0.0, 1.0, (80, 80)).astype(np.float32)
        yy, xx = np.indices(source.shape)
        sigma = 2.2 / 2.355
        source += 18.0 * np.exp(-((xx - 40.2) ** 2 + (yy - 39.8) ** 2) / (2.0 * sigma**2))

        settings = dict(threshold_sigma=4.5, min_pixels=3, expected_fwhm_px=3.0, max_sources=20)
        pixel_detections, _rms = extract_residual_detections(pixel, 0, **settings)
        block_detections, _rms = extract_residual_detections(block, 0, **settings)
        source_detections, _rms = extract_residual_detections(source, 0, **settings)

        self.assertEqual(pixel_detections, [])
        self.assertEqual(block_detections, [])
        self.assertEqual(len(source_detections), 1)
        self.assertGreater(source_detections[0].area_px, 3)

    def test_sensor_fixed_defect_is_removed_after_registration_motion(self):
        defects = {}
        moving_source = {}
        for index in range(6):
            defects[index] = [self._detection(index, 30.0 + index * 1.2, 45.0, 120.0, 80.0)]
            moving_source[index] = [self._detection(index, 30.0 + index * 1.2, 45.0, 120.0 + index * 1.2, 80.0)]

        filtered_defects = suppress_persistent_residuals(defects, persistence_fraction=0.15)
        filtered_source = suppress_persistent_residuals(moving_source, persistence_fraction=0.15)

        self.assertTrue(all(not values for values in filtered_defects.values()))
        self.assertTrue(all(len(values) == 1 for values in filtered_source.values()))

    @staticmethod
    def _detection(frame_index, x, y, sensor_x, sensor_y):
        return Detection(
            frame_index=frame_index,
            x=x,
            y=y,
            ra=None,
            dec=None,
            flux=100.0,
            snr=12.0,
            fwhm_px=2.8,
            eccentricity=0.1,
            area_px=8,
            peak=20.0,
            local_rms=1.0,
            sensor_x=sensor_x,
            sensor_y=sensor_y,
        )


if __name__ == "__main__":
    unittest.main()
