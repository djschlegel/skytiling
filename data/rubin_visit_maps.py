"""Per-band 10-year visit-count maps on a 0.1-deg grid in the NGC and SGC bands from a Rubin opsim visit list
(data/rubin_5.3.5_10yrs.fits reduced to RA, DEC, BAND, PROGRAM, NOTE; a visit covers points within 1.75 deg of its centre;
DD, ToO and twilight visits excluded).  Writes rubin_visits_grid.npy for niji_vs_rubin.py."""
import numpy as np, sys, csv, os
from astropy.io import fits
from scipy.spatial import cKDTree
import os
d = fits.getdata(os.environ.get('RUBIN_VISITS', os.path.expanduser('~/data/djs_tiling/data/subaru/scratch/rubin_5.3.5_10yrs_visits.fits')))
note = np.char.strip(np.array(d['NOTE']))
wfd = ~(np.char.startswith(note, 'DD') | np.char.startswith(note, 'ToO') | np.char.startswith(note, 'twilight'))
print('visits total %d, non-DD/ToO/twilight %d' % (len(d), wfd.sum()))
ra = d['RA'].astype(float); dec = d['DEC'].astype(float); band = np.array(d['BAND'])
def xyz(ra, dec):
    return np.c_[np.cos(np.radians(dec))*np.cos(np.radians(ra)), np.cos(np.radians(dec))*np.sin(np.radians(ra)), np.sin(np.radians(dec))]
R = 1.75; chord = 2*np.sin(np.radians(R)/2)
out = {}
for lab, ralo, rahi, declo, dechi in (('SGC', 320, 410, -10, 10), ('NGC', 122, 232, -6, 8)):
    res = 0.1
    g_ra = (np.arange(ralo, rahi, res) + res/2) % 360; g_dec = np.arange(declo, dechi, res) + res/2
    RA, DEC = np.meshgrid(g_ra, g_dec)
    tree = cKDTree(xyz(RA.ravel(), DEC.ravel()))
    counts = {b: np.zeros(RA.size, np.int32) for b in 'ugrizy'}
    for b in 'ugrizy':
        sel = wfd & (band == b)
        pts = xyz(ra[sel], dec[sel])
        # for each visit, which grid points it covers
        lists = tree.query_ball_point(pts, chord)
        for L in lists:
            counts[b][L] += 1
    out[lab] = dict(RA=RA, DEC=DEC, **{b: counts[b].reshape(RA.shape) for b in 'ugrizy'})
    print(lab, 'done', {b: int(np.median(counts[b][counts[b] > 0])) for b in 'ugrizy'})
np.save('rubin_visits_grid.npy', out, allow_pickle=True)
