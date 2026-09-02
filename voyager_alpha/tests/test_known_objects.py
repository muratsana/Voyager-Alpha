import unittest

import numpy as np

try:
    from astropy.coordinates import SkyCoord  # noqa: F401
except ModuleNotFoundError:
    SkyCoord = None

from voyager_alpha.core.models import Detection


@unittest.skipIf(SkyCoord is None, "astropy dependency is not installed")
class KnownObjectTests(unittest.TestCase):
    def test_match_detection_uses_closest_known_object(self):
        from voyager_alpha.core.known_objects import KnownObject, KnownObjectMatcher

        matcher = KnownObjectMatcher(tolerance_arcsec=30.0)
        detection = Detection(0, 10, 20, 150.0, 22.0, 1000, 10, 3, 0.1)
        objects = [
            KnownObject("far", 151.0, 22.0, None, {}),
            KnownObject("near", 150.001, 22.0, 18.1, {}),
        ]

        match = matcher.match_detection(detection, objects)

        self.assertIsNotNone(match)
        self.assertEqual(match["name"], "near")

    def test_documented_prediction_score_and_status(self):
        from voyager_alpha.core.known_objects import prediction_confidence, prediction_status
        from voyager_alpha.core.models import KnownObjectPrediction

        prediction = KnownObjectPrediction(
            name="test",
            ra=150.0,
            dec=22.0,
            magnitude=17.0,
            local_snr=12.0,
            local_offset_px=1.0,
            expected_trail_px=2.0,
            visible=True,
        )
        prediction.confidence = prediction_confidence(prediction, limiting_magnitude=19.0)

        self.assertGreaterEqual(prediction.confidence, 0.75)
        self.assertEqual(prediction_status(prediction), "high_confidence_match")

    def test_gaia_depth_requires_three_recoveries_per_half_mag_bin(self):
        from astropy.wcs import WCS
        from voyager_alpha.core.known_objects import estimate_gaia_visible_limit

        image = np.random.default_rng(3).normal(100.0, 1.0, (80, 80)).astype(np.float32)
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [40.0, 40.0]
        wcs.wcs.cdelt = [-0.001, 0.001]
        wcs.wcs.crval = [150.0, 22.0]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        rows = []
        for index, magnitude in enumerate((15.1, 15.2, 15.3, 15.6, 15.7)):
            x, y = 20 + index * 9, 35
            yy, xx = np.indices(image.shape)
            image += 30.0 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.2**2))
            ra, dec = wcs.all_pix2world(x, y, 0)
            rows.append({"ra": float(ra), "dec": float(dec), "phot_g_mean_mag": magnitude})

        limit = estimate_gaia_visible_limit(image, wcs.to_header(), catalog_rows=rows)

        self.assertIsNotNone(limit)
        self.assertLess(limit, 15.5)


if __name__ == "__main__":
    unittest.main()
