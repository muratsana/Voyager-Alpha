import numpy as np
from astropy.io import fits
from .fits_io import read_fits_image
import sep

def create_master_background_from_arrays(frames, max_reference_frames=25):
    if not frames:
        raise ValueError("FITS frame listesi boş.")

    if len(frames) > max_reference_frames:
        sample_idx = np.linspace(0, len(frames) - 1, max_reference_frames, dtype=int)
        reference_frames = [frames[i] for i in sample_idx]
    else:
        reference_frames = frames

    first_shape = reference_frames[0].shape
    data_cube = np.zeros((len(reference_frames), *first_shape), dtype=np.float32)
    for index, frame in enumerate(reference_frames):
        if frame.shape != first_shape:
            raise ValueError("FITS boyut uyumsuzluğu.")
        data_cube[index] = np.asarray(frame, dtype=np.float32)
    return np.median(data_cube, axis=0).astype(np.float32)

def create_master_background(fits_files, max_reference_frames=25):
    """
    Verilen FITS dosyalarının medyanını alarak yıldızların sabit kaldığı bir master arka plan üretir.
    Buyuk dizilerde RAM kullanımını sınırlamak için referans kare sayısını sınırlar.
    """
    if not fits_files:
        raise ValueError("FITS dosyası listesi boş.")
        
    first_data = read_fits_image(fits_files[0])
    h, w = first_data.shape
    if len(fits_files) > max_reference_frames:
        sample_idx = np.linspace(0, len(fits_files) - 1, max_reference_frames, dtype=int)
        reference_files = [fits_files[i] for i in sample_idx]
    else:
        reference_files = fits_files

    frames = []
    for file_path in reference_files:
        data = read_fits_image(file_path)
        if data.shape != (h, w):
            raise ValueError(f"FITS boyut uyumsuzluğu: {file_path}")
        frames.append(data)
    return create_master_background_from_arrays(frames, max_reference_frames=max_reference_frames)

def extract_sources(image_data, threshold_sigma=3.0, min_pixels=5):
    """
    SEP kütüphanesini kullanarak fark görüntüsü üzerinden hareketli objeleri (kaynakları) tespit eder.
    SEP native-endian ve C-contiguous float32 veri ister.
    """
    data_sep = np.asarray(image_data, dtype=np.float32)
    if not data_sep.dtype.isnative:
        data_sep = data_sep.byteswap().view(data_sep.dtype.newbyteorder("="))
    if not data_sep.flags['C_CONTIGUOUS']:
        data_sep = np.ascontiguousarray(data_sep)
        
    bkg = sep.Background(data_sep)
    data_sub = data_sep - bkg
    effective_rms = max(float(bkg.globalrms), float(np.nanstd(data_sub)) * 0.25, 1e-6)
    
    objects = sep.extract(data_sub, threshold_sigma, err=effective_rms, minarea=min_pixels)
    return objects, effective_rms
