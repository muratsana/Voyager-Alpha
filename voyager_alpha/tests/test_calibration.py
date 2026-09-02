import unittest

import numpy as np

from voyager_alpha.core.calibration import calibrate_science_frame, calibrate_with_master_dark


class CalibrationTests(unittest.TestCase):
    def test_master_dark_is_scaled_by_exposure(self):
        science = np.full((8, 8), 120.0, dtype=np.float32)
        dark = np.full((8, 8), 10.0, dtype=np.float32)
        calibrated = calibrate_with_master_dark(
            science,
            dark,
            science_exposure=60.0,
            dark_exposure=30.0,
        )
        np.testing.assert_allclose(calibrated, 100.0)

    def test_master_dark_shape_must_match(self):
        with self.assertRaises(ValueError):
            calibrate_with_master_dark(np.zeros((8, 8)), np.zeros((4, 4)))

    def test_bias_dark_and_flat_are_applied_in_order(self):
        bias = np.full((6, 6), 10.0, dtype=np.float32)
        dark = np.full((6, 6), 5.0, dtype=np.float32)
        flat = np.ones((6, 6), dtype=np.float32)
        flat[:, 3:] = 2.0
        signal = np.full((6, 6), 100.0, dtype=np.float32)
        science = signal * (flat / np.median(flat)) + bias + dark * 2.0

        calibrated = calibrate_science_frame(
            science,
            master_bias=bias,
            master_dark=dark,
            master_flat=flat,
            science_exposure=60.0,
            dark_exposure=30.0,
        )

        np.testing.assert_allclose(calibrated, 100.0, atol=1e-5)
