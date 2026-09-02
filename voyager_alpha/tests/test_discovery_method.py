import tempfile
import unittest
from pathlib import Path

import numpy as np

from voyager_alpha.core.discovery_method import DOCUMENTED_DISCOVERY_METHOD
from voyager_alpha.core.pipeline import _temporal_median_from_memmap


class DiscoveryMethodTests(unittest.TestCase):
    def test_documented_defaults_are_centralized(self):
        method = DOCUMENTED_DISCOVERY_METHOD
        self.assertEqual(method.detector_mode, "hybrid")
        self.assertEqual(method.detection_sigma, 5.0)
        self.assertEqual(method.detection_fwhm_px, 3.0)
        self.assertEqual(method.max_residuals_per_frame, 24)
        self.assertEqual(method.edge_margin_px, 6)
        self.assertEqual(method.min_seed_displacement_px, 1.5)
        self.assertEqual(method.match_tolerance_px, 2.8)
        self.assertEqual(method.potential_discovery_rms_px, 0.9)
        self.assertEqual(method.borderline_review_rms_px, 1.8)

    def test_disk_backed_temporal_median_is_exact(self):
        rng = np.random.default_rng(14)
        cube = rng.normal(size=(7, 31, 23)).astype(np.float32)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "cube.dat"
            mapped = np.memmap(path, dtype=np.float32, mode="w+", shape=cube.shape)
            mapped[:] = cube
            mapped.flush()
            actual = _temporal_median_from_memmap(mapped, target_working_mb=1)
            del mapped

        np.testing.assert_allclose(actual, np.median(cube, axis=0), rtol=0.0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
