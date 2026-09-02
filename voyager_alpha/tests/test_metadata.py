import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits

from voyager_alpha.core.metadata import inspect_frame


class FitsMetadataTests(unittest.TestCase):
    def test_reads_camera_sensor_sampling_and_optics_headers(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "camera.fits"
            hdu = fits.PrimaryHDU(np.zeros((40, 60), dtype=np.uint16))
            header = hdu.header
            header["DATE-OBS"] = "2026-09-01T20:00:00"
            header["EXPTIME"] = 180.0
            header["INSTRUME"] = "ZWO ASI2600MM Pro"
            header["DETECTOR"] = "Sony IMX571"
            header["XPIXSZ"] = 3.76
            header["YPIXSZ"] = 3.76
            header["XBINNING"] = 2
            header["YBINNING"] = 2
            header["GAIN"] = 100
            header["OFFSET"] = 50
            header["CCD-TEMP"] = -10.2
            header["SET-TEMP"] = -10.0
            header["FILTER"] = "L"
            header["FOCALLEN"] = 1000.0
            header["APTDIA"] = 200.0
            header["PIXSCALE"] = 1.551
            header["READMODE"] = "High Gain"
            hdu.writeto(path)

            frame = inspect_frame(0, str(path))
            camera = frame.camera

            self.assertEqual(camera.instrument, "ZWO ASI2600MM Pro")
            self.assertEqual(camera.detector, "Sony IMX571")
            self.assertEqual(camera.binning_x, 2)
            self.assertEqual(camera.binning_y, 2)
            self.assertAlmostEqual(camera.pixel_size_x_um, 3.76)
            self.assertAlmostEqual(camera.gain, 100.0)
            self.assertAlmostEqual(camera.offset, 50.0)
            self.assertAlmostEqual(camera.sensor_temperature_c, -10.2)
            self.assertAlmostEqual(camera.image_scale_arcsec_px, 1.551)
            self.assertEqual(camera.image_scale_source, "PIXSCALE")
            self.assertEqual(camera.readout_mode, "High Gain")


if __name__ == "__main__":
    unittest.main()
