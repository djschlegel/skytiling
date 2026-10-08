import numpy as np, os
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
_maps = {}
def _load():
    if not _maps:
        for n, k in ((1, 'ngp'), (-1, 'sgp')):
            _maps[n] = fits.getdata(os.path.join(os.environ.get('SFD_MAPS', os.path.expanduser('~/dustpub/maps')), f'SFD_dust_4096_{k}.fits')).astype(np.float32)
    return _maps
def sfd_ebv(ra, dec):
    """SFD98 E(B-V) at (ra, dec) [deg], bilinear interpolation on the 4096 Lambert maps."""
    m = _load()
    c = SkyCoord(ra=np.asarray(ra, float) * u.deg, dec=np.asarray(dec, float) * u.deg, frame='icrs').galactic
    l = c.l.radian; b = c.b.radian
    out = np.zeros(l.shape, np.float32)
    for n in (1, -1):
        sel = (b >= 0) if n == 1 else (b < 0)
        if not sel.any(): continue
        r = 2048.0 * np.sqrt(1.0 - n * np.sin(b[sel]))
        x = 2048.5 + r * np.cos(l[sel]) - 1.0        # 0-indexed
        y = 2048.5 - n * r * np.sin(l[sel]) - 1.0
        x0 = np.clip(np.floor(x).astype(int), 0, 4094); y0 = np.clip(np.floor(y).astype(int), 0, 4094)
        fx = x - x0; fy = y - y0; img = m[n]
        out[sel] = (img[y0, x0] * (1 - fx) * (1 - fy) + img[y0, x0 + 1] * fx * (1 - fy) + img[y0 + 1, x0] * (1 - fx) * fy + img[y0 + 1, x0 + 1] * fx * fy)
    return out
if __name__ == '__main__':
    # sanity checks: SFD values at known positions
    for name, ra, dec in (('NGP', 192.86, 27.13), ('Coma', 195.0, 28.0), ('M31', 10.68, 41.27), ('Lockman', 161.5, 58.0), ('COSMOS', 150.1, 2.2)):
        print(name, sfd_ebv([ra], [dec]))
