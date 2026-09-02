import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
from astropy.io import fits
from PyQt6.QtWidgets import QApplication

from voyager_alpha.gui.workers import PreviewCacheWorker


class PreviewStfTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_each_blink_frame_receives_its_own_auto_stf(self):
        with tempfile.TemporaryDirectory() as folder:
            files = []
            yy, xx = np.indices((120, 150), dtype=np.float32)
            for index, (background, noise, amplitude) in enumerate(((100.0, 3.0, 240.0), (5000.0, 35.0, 2600.0))):
                rng = np.random.default_rng(500 + index)
                image = rng.normal(background, noise, (120, 150)).astype(np.float32)
                image += amplitude * np.exp(-((xx - 72.0) ** 2 + (yy - 54.0) ** 2) / (2.0 * 1.8**2))
                path = Path(folder) / f"frame_{index}.fits"
                fits.writeto(path, image, overwrite=True)
                files.append(str(path))

            previews = []
            worker = PreviewCacheWorker(files, max_dim=300)
            worker.preview_ready.connect(lambda index, display, shape, thumbnail: previews.append((index, display)))
            worker.run()

            self.assertEqual(len(previews), 2)
            previews.sort(key=lambda item: item[0])
            for _index, display in previews:
                values = np.asarray(display, dtype=np.uint8)
                self.assertGreater(float(np.std(values)), 8.0)
                self.assertGreater(float(np.percentile(values, 99.9) - np.percentile(values, 10.0)), 70.0)

    def test_blink_cache_aligns_shifted_star_field_before_stf(self):
        with tempfile.TemporaryDirectory() as folder:
            yy, xx = np.indices((160, 190), dtype=np.float32)
            stars = [(25, 30), (52, 44), (81, 26), (110, 55), (145, 37), (38, 92), (74, 118), (121, 101), (158, 127)]
            files = []
            for index, (dx, dy, background) in enumerate(((0.0, 0.0, 100.0), (4.2, -3.4, 1700.0))):
                image = np.random.default_rng(900 + index).normal(background, 2.0 + index * 4.0, yy.shape).astype(np.float32)
                for x, y in stars:
                    image += (180.0 + index * 900.0) * np.exp(
                        -((xx - (x + dx)) ** 2 + (yy - (y + dy)) ** 2) / (2.0 * 1.5**2)
                    )
                path = Path(folder) / f"shifted_{index}.fits"
                fits.writeto(path, image, overwrite=True)
                files.append(str(path))

            previews = {}
            alignments = {}
            worker = PreviewCacheWorker(files, max_dim=400)
            worker.preview_ready.connect(lambda index, display, shape, thumbnail: previews.__setitem__(index, display))
            worker.alignment_ready.connect(
                lambda index, method, rms, count: alignments.__setitem__(index, (method, rms, count))
            )
            worker.run()

            offsets = []
            for x, y in stars:
                peaks = []
                for display in (previews[0], previews[1]):
                    patch = display[y - 4 : y + 5, x - 4 : x + 5]
                    peak_y, peak_x = np.unravel_index(int(np.argmax(patch)), patch.shape)
                    peaks.append((x - 4 + peak_x, y - 4 + peak_y))
                offsets.append(float(np.hypot(peaks[1][0] - peaks[0][0], peaks[1][1] - peaks[0][1])))
            self.assertLessEqual(float(np.median(offsets)), 1.0)
            self.assertTrue(alignments[1][0].startswith("star"))
            self.assertGreaterEqual(alignments[1][2], 3)

    def test_blink_cache_corrects_meridian_flip_before_stf(self):
        with tempfile.TemporaryDirectory() as folder:
            yy, xx = np.indices((170, 210), dtype=np.float32)
            reference = np.random.default_rng(77).normal(120.0, 1.2, yy.shape).astype(np.float32)
            stars = ((24, 28), (49, 72), (79, 36), (108, 94), (148, 49), (185, 77), (42, 132), (91, 149), (157, 126))
            for x, y in stars:
                reference += 260.0 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.5**2))
            moving = np.rot90(reference, 2).copy()
            files = []
            for index, image in enumerate((reference, moving)):
                path = Path(folder) / f"flip_{index}.fits"
                fits.writeto(path, image, overwrite=True)
                files.append(str(path))

            previews = {}
            alignments = {}
            worker = PreviewCacheWorker(files, max_dim=400)
            worker.preview_ready.connect(lambda index, display, shape, thumbnail: previews.__setitem__(index, display))
            worker.alignment_ready.connect(
                lambda index, method, rms, count: alignments.__setitem__(index, (method, rms, count))
            )
            worker.run()

            self.assertIn("rot180", alignments[1][0])
            self.assertTrue(alignments[1][0].startswith("star"))
            for x, y in stars:
                first = previews[0][y - 3 : y + 4, x - 3 : x + 4]
                second = previews[1][y - 3 : y + 4, x - 3 : x + 4]
                self.assertLessEqual(
                    float(np.hypot(*(
                        np.asarray(np.unravel_index(int(np.argmax(first)), first.shape), dtype=float)
                        - np.asarray(np.unravel_index(int(np.argmax(second)), second.shape), dtype=float)
                    ))),
                    1.0,
                )


if __name__ == "__main__":
    unittest.main()
