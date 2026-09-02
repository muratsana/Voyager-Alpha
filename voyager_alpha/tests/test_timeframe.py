import unittest

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time

from voyager_alpha.core.timeframe import (
    observatory_location_from_header,
    observing_geometry,
    parse_frame_time,
)


class FrameTimeTests(unittest.TestCase):
    def test_date_avg_is_used_as_midpoint_without_exposure_offset(self):
        header = fits.Header()
        header["DATE-OBS"] = "2026-09-02T20:00:00"
        header["DATE-AVG"] = "2026-09-02T20:00:30"
        header["TIMESYS"] = "UTC"

        result = parse_frame_time(header, 60.0)

        self.assertEqual(result.source, "DATE-AVG")
        self.assertAlmostEqual(result.jd_utc, Time("2026-09-02T20:00:30", scale="utc").jd, places=10)

    def test_date_obs_is_shifted_to_exposure_midpoint(self):
        header = fits.Header()
        header["DATE-OBS"] = "2026-09-02T20:00:00"

        result = parse_frame_time(header, 120.0)

        expected = Time("2026-09-02T20:01:00", scale="utc")
        self.assertAlmostEqual(result.jd_utc, expected.jd, places=10)
        self.assertIn("date_obs_assumed_exposure_start", result.flags)

    def test_barycentric_time_matches_astropy_reference_better_than_one_ms(self):
        header = fits.Header()
        header["SITELAT"] = 39.9334
        header["SITELONG"] = 32.8597
        header["SITEELEV"] = 900.0
        location = observatory_location_from_header(header)
        target = SkyCoord(ra=120.5 * u.deg, dec=22.3 * u.deg)
        midpoint = Time("2026-09-02T20:00:30", scale="utc")

        result = observing_geometry(midpoint, target, location)

        reference_time = Time(midpoint.jd, format="jd", scale="utc", location=EarthLocation.from_geodetic(32.8597, 39.9334, 900.0))
        expected = (reference_time.tdb + reference_time.light_travel_time(target, kind="barycentric")).jd
        error_ms = abs(result.bjd_tdb - expected) * 86400.0 * 1000.0
        self.assertLess(error_ms, 1.0)
        self.assertIsNotNone(result.altitude_deg)
        self.assertIsNotNone(result.sun_altitude_deg)


if __name__ == "__main__":
    unittest.main()
