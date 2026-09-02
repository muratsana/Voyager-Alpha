import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from voyager_alpha.core.exoplanet_catalog import ExoplanetCatalog, predict_nearest_transit
from voyager_alpha.core.exoplanet_worker import SequencePlateSolveWorker
from voyager_alpha.core.models import RegistrationSolution
from voyager_alpha.core.plate_solver import PlateSolveResult
from voyager_alpha.core.wcs_cache import WcsSolutionCache, pixel_to_sky, propagate_wcs_header


class ExoplanetCatalogTests(unittest.TestCase):
    def test_updates_all_supported_sources_and_searches_by_position(self):
        payloads = {
            "nea_confirmed": (
                "pl_name,hostname,ra,dec,pl_orbper,pl_tranmid,pl_trandur,pl_trandep,disc_facility\n"
                "Test b,Test,120.0,22.0,2.0,2460000.0,2.4,1.2,Transit Survey\n"
            ),
            "tess_toi": (
                "toidisplay,tid,toi,ra,dec,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,pl_trandep,st_tmag,rowupdate\n"
                "TOI-100.01,123,100.01,120.0002,22.0001,PC,3.0,2460001.0,3.1,1500,10.2,2026-01-01\n"
            ),
            "kepler_koi": (
                "kepoi_name,kepler_name,kepid,ra,dec,koi_disposition,koi_period,koi_time0bk,koi_duration,koi_depth,koi_kepmag,koi_vet_date\n"
                "K00001.01,,42,120.0003,22.0,FALSE POSITIVE,4.0,100.0,2.0,900,14.2,2025-01-01\n"
            ),
            "k2_candidates": (
                "pl_name,epic_candname,k2_name,hostname,ra,dec,disposition,pl_orbper,pl_tranmid,pl_trandur,pl_trandep,sy_kepmag,rowupdate\n"
                "EPIC 1.01,EPIC 1.01,,EPIC 1,120.0004,22.0,CANDIDATE,5.0,2460002.0,4.0,0.5,12.1,2024-01-01\n"
            ),
        }
        with tempfile.TemporaryDirectory() as folder:
            catalog = ExoplanetCatalog(Path(folder) / "catalog.sqlite3")
            counts = catalog.update_all(downloader=lambda source: payloads[source.key])
            self.assertEqual(sum(counts.values()), 4)
            self.assertEqual(catalog.total_records(), 4)
            matches = catalog.cone_search(120.0, 22.0, radius_arcmin=1.0)
            self.assertEqual({item.disposition for item in matches}, {"confirmed", "candidate", "false_positive"})
            koi = next(item for item in matches if item.source_key == "kepler_koi")
            self.assertAlmostEqual(koi.epoch_bjd, 2454933.0)
            confirmed = next(item for item in matches if item.source_key == "nea_confirmed")
            prediction = predict_nearest_transit(confirmed, 2459999.95, 2460000.05)
            self.assertIsNotNone(prediction)
            self.assertTrue(prediction.overlaps_observation)

    def test_wcs_cache_and_affine_propagation_preserve_sky_position(self):
        reference = WCS(naxis=2)
        reference.wcs.crpix = [50.0, 40.0]
        reference.wcs.cdelt = [-1.2 / 3600.0, 1.2 / 3600.0]
        reference.wcs.crval = [150.0, 25.0]
        reference.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        reference_header = reference.to_header()
        registration = RegistrationSolution(
            matrix_xy=np.asarray([[1.001, 0.002], [-0.001, 0.999]], dtype=float),
            offset_xy=np.asarray([4.5, -2.0], dtype=float),
            rms_px=0.2,
            matched_stars=30,
            method="star-affine",
        )
        moving_header = propagate_wcs_header(reference_header, registration)
        moving_xy = np.asarray([32.0, 61.0])
        reference_xy = registration.matrix_xy @ moving_xy + registration.offset_xy
        expected = pixel_to_sky(reference_header, *reference_xy)
        actual = pixel_to_sky(moving_header, *moving_xy)
        self.assertAlmostEqual(expected[0], actual[0], places=8)
        self.assertAlmostEqual(expected[1], actual[1], places=8)

        with tempfile.TemporaryDirectory() as folder:
            fits_path = Path(folder) / "frame.fits"
            fits.writeto(fits_path, np.zeros((80, 100), dtype=np.float32), overwrite=True)
            cache = WcsSolutionCache(Path(folder) / "wcs")
            cache.store(str(fits_path), moving_header)
            loaded, method = cache.resolve(str(fits_path))
            self.assertEqual(method, "cache")
            self.assertIsNotNone(loaded)

    def test_sequence_plate_solver_solves_reference_and_propagates_other_frames(self):
        with tempfile.TemporaryDirectory() as folder:
            files = []
            yy, xx = np.indices((140, 160), dtype=np.float32)
            stars = [(20 + (index % 5) * 27, 20 + (index // 5) * 32) for index in range(15)]
            for frame_index in range(3):
                image = np.full((140, 160), 100.0, dtype=np.float32)
                for star_index, (x, y) in enumerate(stars):
                    image += (800 + star_index * 25) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 1.5**2))
                path = Path(folder) / f"frame_{frame_index}.fits"
                fits.writeto(path, image, overwrite=True)
                files.append(str(path))

            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [80.0, 70.0]
            wcs.wcs.cdelt = [-1.4 / 3600.0, 1.4 / 3600.0]
            wcs.wcs.crval = [180.0, 30.0]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

            class FakeSolver:
                calls = 0

                @staticmethod
                def is_available():
                    return True

                def solve(self, _path):
                    self.calls += 1
                    return PlateSolveResult(True, header=wcs.to_header())

            solver = FakeSolver()
            cache = WcsSolutionCache(Path(folder) / "wcs-cache")
            results = []
            worker = SequencePlateSolveWorker(files, cache=cache, solver=solver)
            worker.result_ready.connect(results.append)
            worker.run()

            self.assertEqual(len(results), 1)
            self.assertTrue(results[0].complete)
            self.assertEqual(solver.calls, 1)
            self.assertTrue(all(cache.resolve(path)[0] is not None for path in files))
            self.assertGreaterEqual(sum(method.startswith("propagated") for method in results[0].methods.values()), 2)


if __name__ == "__main__":
    unittest.main()
