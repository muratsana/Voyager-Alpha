import unittest

from voyager_alpha.core.models import Detection, FrameRecord, Tracklet
from voyager_alpha.core.reporting import render_tracklet_html_report


class ReportingTests(unittest.TestCase):
    def test_html_report_contains_tracklet_and_quality_summary(self):
        frames = [
            FrameRecord(0, "a.fits", "2025-01-01T00:00:00", 2460676.5, 60, (50, 50), True),
            FrameRecord(1, "b.fits", "2025-01-01T00:01:00", 2460676.5007, 60, (50, 50), True),
            FrameRecord(2, "c.fits", "2025-01-01T00:02:00", 2460676.5014, 60, (50, 50), True),
        ]
        detections = [
            Detection(0, 10, 20, None, None, 1000, 9, 3, 0.1),
            Detection(1, 13, 22, None, None, 1000, 10, 3, 0.1),
            Detection(2, 16, 24, None, None, 1000, 11, 3, 0.1),
        ]
        tracklet = Tracklet("AH-00001", detections, 3, 3.6, 4.3, 120.0, 0.2, 10.0, 0.8)

        html = render_tracklet_html_report(frames, [tracklet])

        self.assertIn("Voyager Alpha Asteroid Report", html)
        self.assertIn("AH-00001", html)
        self.assertIn("WCS OK", html)


if __name__ == "__main__":
    unittest.main()
