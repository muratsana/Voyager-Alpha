import unittest
import tempfile
from pathlib import Path

from voyager_alpha.core.plate_solver import AstapPlateSolver, _looks_solved, _read_astap_status
from voyager_alpha.core.models import FrameRecord
from voyager_alpha.core.pipeline import _frame_fov_deg


class PlateSolverTests(unittest.TestCase):
    def test_reports_unavailable_without_executable(self):
        solver = AstapPlateSolver(executable="Z:/missing/astap.exe")

        self.assertFalse(solver.is_available())
        self.assertFalse(solver.solve("missing.fits").success)

    def test_output_marker_detection(self):
        self.assertTrue(_looks_solved("Solution found, WCS written"))
        self.assertFalse(_looks_solved("No solution found"))

    def test_reads_astap_machine_status(self):
        with tempfile.TemporaryDirectory() as folder:
            working = Path(folder) / "frame.fits"
            working.with_suffix(".ini").write_text(
                "PLTSOLVD=T\nWARNING=Field of view estimate refined\n",
                encoding="utf-8",
            )

            status, warning, log_path = _read_astap_status(working, Path(folder))

            self.assertEqual(status, "T")
            self.assertIn("refined", warning)
            self.assertIsNone(log_path)

    def test_astap_fov_uses_vertical_sensor_dimension(self):
        frame = FrameRecord(0, "frame.fits", None, None, 60.0, (2000, 3000), False)
        frame.camera.image_scale_arcsec_px = 1.8

        self.assertAlmostEqual(_frame_fov_deg(frame), 1.0)


if __name__ == "__main__":
    unittest.main()
