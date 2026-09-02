import unittest

from voyager_alpha.core.ades import render_ades_psv_draft
from voyager_alpha.core.models import Detection, FrameRecord, Tracklet


class AdesDraftTests(unittest.TestCase):
    def test_renders_draft_warning_and_tracklet_rows(self):
        frame = FrameRecord(0, "a.fits", "2025-01-01T00:00:00", 2460676.5, 60, (50, 50), True)
        detection = Detection(0, 10, 20, 150.1, 22.2, 1000, 9, 3, 0.1)
        tracklet = Tracklet("AH-00001", [detection], 1, 1.0, 1.2, 90.0, 0.1, 9.0, 0.5)

        psv = render_ades_psv_draft([tracklet], [frame], observatory_code="XXX")

        self.assertIn("trkSub", psv)
        self.assertIn("AH-00001", psv)
        self.assertIn("DRAFT_NOT_MPC_READY", psv)
        self.assertIn("XXX", psv)


if __name__ == "__main__":
    unittest.main()

