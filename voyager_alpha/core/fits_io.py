from astropy.io import fits


def read_fits_image(file_path: str, *, header: bool = False):
    """Read scaled FITS image data safely for BSCALE/BZERO/BLANK files."""
    data, hdr = fits.getdata(file_path, header=True, memmap=False)
    if header:
        return data, hdr
    return data


def read_fits_header(file_path: str):
    return fits.getheader(file_path, 0)
