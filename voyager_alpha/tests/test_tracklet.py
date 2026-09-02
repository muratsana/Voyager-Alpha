import unittest

from voyager_alpha.core.models import Detection, FrameRecord
from voyager_alpha.core.detection import suppress_persistent_residuals
from voyager_alpha.core.tracklet import link_tracklets


def frame(index):
    return FrameRecord(
        index=index,
        file_path=f"f{index}.fits",
        date_obs="2025-01-01T00:00:00",
        midpoint_jd=2460676.5 + index / 1440.0,
        exposure_seconds=60.0,
        shape=(100, 100),
        has_wcs=True,
    )


def detection(frame_index, x, y, snr=10.0):
    return Detection(
        frame_index=frame_index,
        x=x,
        y=y,
        ra=None,
        dec=None,
        flux=1000.0,
        snr=snr,
        fwhm_px=3.0,
        eccentricity=0.1,
    )


class TrackletTests(unittest.TestCase):
    def test_persistent_filter_removes_fixed_residual_but_keeps_slow_motion(self):
        fixed = {
            index: [detection(index, 25.0 + (index % 2) * 0.05, 30.0)]
            for index in range(20)
        }
        slow = {
            index: [detection(index, 10.0 + index * 0.35, 40.0 + index * 0.08)]
            for index in range(20)
        }

        fixed_filtered = suppress_persistent_residuals(fixed, persistence_fraction=0.15)
        slow_filtered = suppress_persistent_residuals(slow, persistence_fraction=0.15)

        self.assertEqual(sum(map(len, fixed_filtered.values())), 0)
        self.assertEqual(sum(map(len, slow_filtered.values())), 20)

    def test_links_linear_motion_across_frames(self):
        frames = [frame(i) for i in range(4)]
        detections = {i: [detection(i, 10 + i * 3, 20 + i * 2)] for i in range(4)}

        tracklets = link_tracklets(detections, frames, pixel_scale_arcsec=1.0)

        self.assertEqual(len(tracklets), 1)
        self.assertEqual(tracklets[0].frames_detected, 4)
        self.assertGreater(tracklets[0].confidence, 0.8)

    def test_rejects_stationary_hot_pixel(self):
        frames = [frame(i) for i in range(4)]
        detections = {i: [detection(i, 25, 25)] for i in range(4)}

        tracklets = link_tracklets(detections, frames, pixel_scale_arcsec=1.0)

        self.assertEqual(tracklets, [])

    def test_rejects_low_snr_tracklet(self):
        frames = [frame(i) for i in range(4)]
        detections = {i: [detection(i, 10 + i * 3, 20 + i * 2, snr=3.5)] for i in range(4)}

        tracklets = link_tracklets(detections, frames, pixel_scale_arcsec=1.0)

        self.assertEqual(tracklets, [])

    def test_links_constant_velocity_with_irregular_cadence(self):
        minute_offsets = [0.0, 0.7, 2.4, 4.1, 7.0]
        frames = []
        detections = {}
        for index, minute in enumerate(minute_offsets):
            record = frame(index)
            record.midpoint_jd = 2460676.5 + minute / 1440.0
            frames.append(record)
            detections[index] = [detection(index, 18.0 + 2.5 * minute, 22.0 - 1.2 * minute)]

        tracklets = link_tracklets(detections, frames, min_frames=4, pixel_scale_arcsec=1.0)

        self.assertEqual(len(tracklets), 1)
        self.assertEqual(tracklets[0].frames_detected, 5)
        self.assertAlmostEqual(tracklets[0].velocity_x_px_min, 2.5, places=1)
        self.assertAlmostEqual(tracklets[0].velocity_y_px_min, -1.2, places=1)

    def test_rejects_sparse_random_chain_across_long_span(self):
        frames = [frame(index) for index in range(12)]
        detections = {index: [] for index in range(12)}
        detections[0] = [detection(0, 10, 10)]
        detections[1] = [detection(1, 12, 11)]
        detections[10] = [detection(10, 30, 20)]
        detections[11] = [detection(11, 32, 21)]

        tracklets = link_tracklets(
            detections,
            frames,
            min_frames=4,
            min_track_occupancy=0.6,
            max_missing_gap_frames=2,
        )

        self.assertEqual(tracklets, [])

    def test_rejects_artifact_dominated_linear_chain(self):
        frames = [frame(index) for index in range(5)]
        detections = {}
        for index in range(5):
            item = detection(index, 15 + index * 2, 20 + index)
            item.flags.append("star_subtraction_dipole")
            detections[index] = [item]

        tracklets = link_tracklets(detections, frames, max_artifact_fraction=0.34)

        self.assertEqual(tracklets, [])


if __name__ == "__main__":
    unittest.main()
