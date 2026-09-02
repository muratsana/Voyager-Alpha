import unittest

import numpy as np
from scipy.ndimage import shift

from voyager_alpha.core.registration import estimate_integer_shift, estimate_subpixel_shift, register_frame, shift_image_integer


class RegistrationTests(unittest.TestCase):
    def test_estimates_shift_needed_to_align_moving_frame(self):
        reference = np.zeros((64, 64), dtype=np.float32)
        reference[20, 30] = 100.0
        reference[40, 10] = 60.0

        moving = shift_image_integer(reference, 4, -3)
        dy, dx, peak = estimate_integer_shift(reference, moving)

        self.assertEqual((dy, dx), (-4, 3))
        self.assertGreater(peak, 0.9)

    def test_estimates_subpixel_shift(self):
        yy, xx = np.indices((96, 96), dtype=np.float32)
        reference = np.exp(-((xx - 33.4) ** 2 + (yy - 47.2) ** 2) / (2 * 2.1**2))
        reference += 0.7 * np.exp(-((xx - 68.3) ** 2 + (yy - 25.8) ** 2) / (2 * 1.7**2))
        moving = shift(reference, shift=(2.35, -1.65), order=3, mode="constant", cval=0.0)

        dy, dx, peak = estimate_subpixel_shift(reference, moving)

        self.assertAlmostEqual(dy, -2.35, delta=0.45)
        self.assertAlmostEqual(dx, 1.65, delta=0.45)
        self.assertGreater(peak, 0.1)

    def test_uses_star_translation_when_three_to_seven_stars_match(self):
        yy, xx = np.indices((120, 140), dtype=np.float32)
        reference = np.zeros_like(xx)
        for x, y in ((22, 20), (48, 76), (73, 35), (101, 91), (119, 52)):
            reference += 100.0 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.5**2))
        moving = shift(reference, shift=(3.25, -2.5), order=3, mode="constant", cval=0.0)

        aligned, solution = register_frame(reference, moving)
        dy, dx, _peak = estimate_subpixel_shift(reference, aligned)

        self.assertEqual(solution.method, "star-translation")
        self.assertGreaterEqual(solution.matched_stars, 3)
        self.assertAlmostEqual(dx, 0.0, delta=0.45)
        self.assertAlmostEqual(dy, 0.0, delta=0.45)

    def test_corrects_meridian_flip_before_star_alignment(self):
        yy, xx = np.indices((180, 220), dtype=np.float32)
        reference = np.random.default_rng(42).normal(100.0, 0.8, yy.shape).astype(np.float32)
        stars = (
            (21.4, 24.8, 180.0),
            (48.2, 62.7, 260.0),
            (78.8, 31.5, 145.0),
            (109.1, 88.6, 310.0),
            (151.7, 43.2, 220.0),
            (190.3, 71.4, 170.0),
            (37.6, 126.2, 280.0),
            (86.5, 147.1, 195.0),
            (137.9, 118.4, 240.0),
            (181.2, 153.3, 330.0),
        )
        for x, y, amplitude in stars:
            reference += amplitude * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.45**2))
        moving = shift(np.rot90(reference, 2), shift=(3.35, -4.6), order=3, mode="constant", cval=100.0)

        aligned, solution = register_frame(reference, moving)
        dy, dx, _peak = estimate_subpixel_shift(reference, aligned)

        self.assertIn("rot180", solution.method)
        self.assertTrue(solution.method.startswith("star"))
        self.assertGreaterEqual(solution.matched_stars, 8)
        self.assertAlmostEqual(dx, 0.0, delta=0.25)
        self.assertAlmostEqual(dy, 0.0, delta=0.25)

    def test_normal_affine_alignment_leaves_subpixel_star_residual(self):
        yy, xx = np.indices((180, 220), dtype=np.float32)
        reference = np.full(yy.shape, 50.0, dtype=np.float32)
        for x, y in ((24, 26), (51, 65), (82, 34), (111, 91), (154, 46), (192, 73), (39, 129), (88, 149), (140, 121)):
            reference += 220.0 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.4**2))
        moving = shift(reference, shift=(-2.65, 3.4), order=3, mode="constant", cval=50.0)

        aligned, solution = register_frame(reference, moving)
        dy, dx, _peak = estimate_subpixel_shift(reference, aligned)

        self.assertNotIn("rot180", solution.method)
        self.assertTrue(solution.method.startswith("star"))
        self.assertAlmostEqual(dx, 0.0, delta=0.2)
        self.assertAlmostEqual(dy, 0.0, delta=0.2)


if __name__ == "__main__":
    unittest.main()
