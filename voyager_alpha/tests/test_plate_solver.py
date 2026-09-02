import unittest

from voyager_alpha.core.plate_solver import AstapPlateSolver, _looks_solved


class PlateSolverTests(unittest.TestCase):
    def test_reports_unavailable_without_executable(self):
        solver = AstapPlateSolver(executable="Z:/missing/astap.exe")

        self.assertFalse(solver.is_available())
        self.assertFalse(solver.solve("missing.fits").success)

    def test_output_marker_detection(self):
        self.assertTrue(_looks_solved("Solution found, WCS written"))
        self.assertFalse(_looks_solved("No solution found"))


if __name__ == "__main__":
    unittest.main()

