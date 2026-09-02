import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    from astropy.io import fits
    from astropy.wcs import WCS
    import sep  # noqa: F401
except ModuleNotFoundError:
    fits = None
    WCS = None

from voyager_alpha.core.models import Detection
from voyager_alpha.core.tracklet import link_tracklets


@unittest.skipIf(fits is None or WCS is None, "astropy/sep dependencies are not installed")
class SyntheticPipelineTests(unittest.TestCase):
    def test_synthetic_moving_object_becomes_tracklet(self):
        from voyager_alpha.core.image_math import create_master_background_from_arrays, extract_sources
        from voyager_alpha.core.metadata import inspect_sequence

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            files = []
            frames = []
            for index in range(5):
                image = self._synthetic_frame(index)
                frames.append(image)
                header = self._wcs_header(image.shape, index)
                file_path = tmp_path / f"frame_{index:03d}.fits"
                fits.writeto(file_path, image.astype(np.float32), header, overwrite=True)
                files.append(str(file_path))

            frame_records = inspect_sequence(files)
            master = create_master_background_from_arrays(frames)
            detections_by_frame = {}
            for index, image in enumerate(frames):
                diff = np.clip(image - master, 0, None)
                objects, global_rms = extract_sources(diff, threshold_sigma=3.0, min_pixels=2)
                detections_by_frame[index] = [
                    Detection(
                        frame_index=index,
                        x=float(obj["x"]),
                        y=float(obj["y"]),
                        ra=None,
                        dec=None,
                        flux=float(obj["flux"]),
                        snr=float(obj["flux"] / max(global_rms * np.sqrt(max(float(obj["npix"]), 1.0)), 1e-6)),
                        fwhm_px=3.0,
                        eccentricity=0.1,
                    )
                    for obj in objects
                ]

            tracklets = link_tracklets(detections_by_frame, frame_records, min_frames=3, pixel_scale_arcsec=1.2)

            self.assertGreaterEqual(len(tracklets), 1)
            self.assertGreaterEqual(tracklets[0].frames_detected, 3)

    def test_streaming_worker_completes_sequence(self):
        from voyager_alpha.core.pipeline import AsteroidWorker

        with tempfile.TemporaryDirectory() as tmp:
            files = []
            for index in range(5):
                image = self._synthetic_frame(index)
                path = Path(tmp) / f"stream_{index:03d}.fits"
                fits.writeto(path, image, self._wcs_header(image.shape, index), overwrite=True)
                files.append(str(path))

            sequences = []
            frame_results = []
            results = []
            worker = AsteroidWorker(
                files,
                sigma=3.0,
                min_pix=2,
                min_tracklet_frames=3,
                match_known_objects=False,
                estimate_gaia_depth=False,
            )
            worker.sequence_ready.connect(sequences.append)
            worker.frame_done.connect(lambda index, result: frame_results.append((index, result)))
            worker.result_ready.connect(results.append)
            worker.run()

            self.assertTrue(sequences)
            self.assertEqual(len(frame_results), 5)
            self.assertEqual(len(results), 1)
            self.assertGreaterEqual(len(results[0].tracklets), 1)

    def _synthetic_frame(self, index):
        rng = np.random.default_rng(100 + index)
        image = rng.normal(100.0, 1.5, size=(96, 96)).astype(np.float32)
        static_stars = [(15, 20), (30, 70), (70, 45), (80, 80)]
        for y, x in static_stars:
            self._add_gaussian(image, x, y, amplitude=250.0)
        self._add_gaussian(image, 25 + index * 3, 35 + index * 2, amplitude=180.0)
        return image

    def _add_gaussian(self, image, x, y, amplitude):
        yy, xx = np.indices(image.shape)
        image += amplitude * np.exp(-(((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 1.2**2)))

    def _wcs_header(self, shape, index):
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [shape[1] / 2, shape[0] / 2]
        wcs.wcs.cdelt = np.array([-0.00033, 0.00033])
        wcs.wcs.crval = [150.0, 22.0]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        header = wcs.to_header()
        header["DATE-OBS"] = f"2025-01-01T00:0{index}:00"
        header["EXPTIME"] = 60.0
        return header


if __name__ == "__main__":
    unittest.main()
