from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import warnings

def pixel_to_radec(header, x, y):
    """WCS bilgisi içeren FITS başlığını kullanarak pikseli RA/DEC'e çevirir."""
    try:
        wcs_info = WCS(header)
        if wcs_info.has_celestial:
            ra, dec = wcs_info.all_pix2world(x, y, 0)
            return float(ra), float(dec)
    except Exception:
        pass
    return None, None

def estimate_pixel_scale_arcsec(header):
    """FITS WCS bilgisinden yaklaşık piksel ölçeğini arcsec/pixel olarak döndürür."""
    try:
        wcs_info = WCS(header)
        if not wcs_info.has_celestial:
            return None
        scales_deg = proj_plane_pixel_scales(wcs_info.celestial)
        return float(abs(scales_deg.mean()) * 3600.0)
    except Exception:
        return None

def check_known_asteroid(ra, dec, date_obs_str, search_radius=60, tolerance=15):
    """SkyBoT (MPC) veritabanında verilen koordinat ve zamanda asteroit olup olmadığını sorgular."""
    if not date_obs_str or ra is None or dec is None:
        return False, None
        
    try:
        from astropy.time import Time
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        from astroquery.imcce import Skybot

        epoch = Time(date_obs_str, format='fits')
        target_coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            results = Skybot.cone_search(target_coord, rad=search_radius*u.arcsec, epoch=epoch)
            
        if results is None or len(results) == 0:
            return False, None
            
        min_sep = tolerance
        best_match = None
        
        for row in results:
            db_coord = SkyCoord(ra=row['RA']*u.degree, dec=row['DEC']*u.degree, frame='icrs')
            sep_arcsec = target_coord.separation(db_coord).arcsec
            
            if sep_arcsec <= min_sep:
                min_sep = sep_arcsec
                best_match = {'name': str(row['Name']), 'mag': float(row['V']), 'sep': sep_arcsec}
                
        if best_match:
            return True, best_match
            
    except Exception:
        pass
        
    return False, None
