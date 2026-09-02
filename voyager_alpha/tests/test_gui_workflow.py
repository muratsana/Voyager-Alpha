import os
import tempfile
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
from astropy.io import fits
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QFrame

from voyager_alpha.gui.main_window import MainWindow
from voyager_alpha.main import APP_USER_MODEL_ID, application_icon_path


class GuiWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_two_module_shell_and_scientific_workspaces(self):
        with tempfile.TemporaryDirectory() as folder:
            rng = np.random.default_rng(42)
            for index in range(5):
                data = rng.normal(1000, 12, (96, 128)).astype(np.float32)
                self._add_star(data, 55, 40, 900 - index * 5)
                self._add_star(data, 72, 48, 800)
                hdu = fits.PrimaryHDU(data)
                hdu.header["DATE-OBS"] = f"2026-09-01T20:0{index}:00"
                hdu.header["EXPTIME"] = 60.0
                hdu.header["INSTRUME"] = "Test Astro Camera"
                hdu.header["DETECTOR"] = "IMX Test Sensor"
                hdu.header["XPIXSZ"] = 3.76
                hdu.header["YPIXSZ"] = 3.76
                hdu.header["XBINNING"] = 1
                hdu.header["YBINNING"] = 1
                hdu.header["GAIN"] = 100
                hdu.header["CCD-TEMP"] = -10.0
                hdu.writeto(os.path.join(folder, f"frame_{index + 1:03d}.fits"))

            window = MainWindow()
            try:
                window.resize(1920, 1080)
                window.show()
                self.app.processEvents()
                self.assertEqual(window.windowTitle(), "Voyager Alpha")
                self.assertFalse(window.windowIcon().isNull())
                self.assertEqual(APP_USER_MODEL_ID, "tr.com.astrohub.voyageralpha")
                self.assertEqual(application_icon_path().name, "voyager-alpha.ico")
                self.assertTrue(application_icon_path().is_file())
                self.assertEqual(window.lbl_brand_name.text(), "Voyager Alpha")
                self.assertIsNotNone(window.findChild(QFrame, "logoPlate"))
                self.assertEqual(window.stack.count(), 2)
                self.assertEqual(window.btn_asteroid.text(), "Asteroid Hunter")
                self.assertEqual(window.btn_exoplanet.text(), "Exoplanet Inspection")

                asteroid = window.asteroid_workspace
                asteroid.scan_folder(folder)
                self.app.processEvents()
                self.assertEqual(len(asteroid.fits_files), 5)
                self.assertEqual(len(asteroid.filmstrip.paths), 5)
                camera_rows = asteroid.sequence_table.findItems("Camera", Qt.MatchFlag.MatchExactly)
                self.assertEqual(len(camera_rows), 1)
                self.assertEqual(asteroid.sequence_table.item(camera_rows[0].row(), 1).text(), "Test Astro Camera")
                self.assertEqual(
                    asteroid.sequence_table.verticalScrollBarPolicy(),
                    Qt.ScrollBarPolicy.ScrollBarAlwaysOff,
                )
                self.assertEqual(asteroid.filmstrip.height(), 108)
                self.assertEqual(
                    asteroid.tracklet_table.horizontalScrollBarPolicy(),
                    Qt.ScrollBarPolicy.ScrollBarAlwaysOff,
                )
                self.assertGreaterEqual(asteroid.lbl_file_status.width(), 300)
                self.assertLessEqual(asteroid.pipeline_progress.width(), 720)
                self.assertTrue(asteroid.btn_optimize_detection.isEnabled())
                asteroid.combo_blink_speed.setCurrentText("10 fps")
                self.assertEqual(asteroid.current_blink_interval(), 100)
                window.resize(1280, 760)
                self.app.processEvents()
                self.assertLess(asteroid.filmstrip.page_size, asteroid.filmstrip.max_page_size)
                self.assertGreaterEqual(asteroid.filmstrip.page_size, 3)

                window.btn_exoplanet.click()
                self.assertEqual(window.stack.currentWidget(), window.exoplanet_workspace)
                exoplanet = window.exoplanet_workspace
                self.assertEqual(exoplanet.right_tabs.count(), 4)
                self.assertEqual(exoplanet.btn_auto_references.text(), "Otomatik Referans")
                self.assertEqual(exoplanet.spin_linearity.value(), 0.0)
                self.assertIn("ÖN KONTROL", exoplanet.lbl_quality_verdict.text())
                self.assertEqual(exoplanet.catalog_sources_table.rowCount(), 4)
                self.assertTrue(exoplanet.btn_plate_solve.isEnabled() is False)
                exoplanet.scan_folder(folder)
                self.app.processEvents()
                self.assertEqual(
                    exoplanet.validation.verticalScrollBarPolicy(),
                    Qt.ScrollBarPolicy.ScrollBarAlwaysOff,
                )
                self.assertGreaterEqual(exoplanet.lbl_file_status.width(), 300)
                self.assertTrue(exoplanet.btn_plate_solve.isEnabled())
                exoplanet.set_selection_mode("target")
                exoplanet.image_clicked(55.0, 40.0)
                exoplanet.set_selection_mode("comparison")
                exoplanet.image_clicked(72.0, 48.0)
                self.assertEqual(exoplanet.target_xy, (55.0, 40.0))
                self.assertEqual(len(exoplanet.comparison_xy), 1)
                self.assertTrue(exoplanet.btn_start.isEnabled())
            finally:
                window.close()

    @staticmethod
    def _add_star(image, x, y, amplitude):
        yy, xx = np.indices(image.shape)
        image += amplitude * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 1.7**2))


if __name__ == "__main__":
    unittest.main()
