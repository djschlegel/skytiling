"""Evaluate stripe tilings on independent randoms: inside/interior stats, spill, Dec edge profiles."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/repo')
from skytiling.optimize_tiling_maps import *
from skytiling.optimize_tiling_maps import _tile_coverage
from astropy.io import fits
from scipy.spatial import cKDTree
cam = ThroughputCamera('/home/claude/repo/data/subaru3/hsc_mbq1_throughput.fits')
fp = parse_rects('/home/claude/repo/data/subaru3/niji_stripes_2p8_hsc.csv')
rng = np.random.default_rng(11)
rand = fp.random_points(6_000_000, rng, extra=1.2)          # footprint + 1.2 deg so spill is measured
tree = cKDTree(rand); inside = fp.contains(rand); interior = fp.interior(rand, 0.75)
ra, dec = xyz_to_radec(rand)
rot = [np.radians(r) for r in (0, 90, 180, 270)]; chord = 2*np.sin(np.radians(cam.max_radius)/2)
def coverage(tt):
    tiles = radec_to_xyz(tt['ra'], tt['dec']); theta = np.radians(tt['rot'])
    c = np.zeros((len(rand), cam.nf), np.float32)
    for i in range(len(tiles)):
        k = np.array(tree.query_ball_point(tiles[i], chord), dtype=np.int64)
        if len(k): c[k] += _tile_coverage(cam, rot, tiles[i], rand[k], theta[i])
    return c
def stats(c, sel):
    cc = c[sel]; m = cc.mean(0); s = cc.std(0)
    return m, s, 100*(cc < 1.5).mean(0), 100*(cc < 0.5).mean(0)
runs = dict(sys.argv[i].split('=') for i in range(1, len(sys.argv)))
import os
out = {}
for name, path in runs.items():
    tt = fits.getdata(path); c = coverage(tt)
    tot = c.sum(0); spill = 100*c[~inside].sum(0)/tot
    m, s, h15, h05 = stats(c, inside); mi, si, hi15, hi05 = stats(c, interior)
    print(f'{name:12s} spill outside {np.round(spill,1)}%')
    print(f'   inside   mean {np.round(m,2)} rms/mean {np.round(s/m,3)} <1.5 {np.round(h15,1)}% <0.5 {np.round(h05,2)}%')
    print(f'   interior mean {np.round(mi,2)} rms/mean {np.round(si/mi,3)} <1.5 {np.round(hi15,1)}% <0.5 {np.round(hi05,2)}%', flush=True)
    # Dec profile across the NGC band (RA 140-215) and the SGC band (RA 340-30), 413 and mean of 4
    prof = {}
    for lab, rsel, dlo, dhi in (('NGC', (ra > 140) & (ra < 215), -3.5, 4.9), ('SGC', ((ra > 340) | (ra < 28)), -2.3, 6.1)):
        bins = np.arange(dlo, dhi + 0.001, 0.1); idx = np.digitize(dec, bins) - 1
        p = np.full((len(bins)-1, cam.nf), np.nan)
        for b in range(len(bins)-1):
            ss = rsel & (idx == b)
            if ss.sum() > 50: p[b] = c[ss].mean(0)
        prof[lab] = (0.5*(bins[1:]+bins[:-1]), p)
    out[name] = dict(spill=spill, inside=(m, s, h15, h05), interior=(mi, si, hi15, hi05), prof=prof, cov_in=c[inside, 0][::20].copy())
np.save(sys.argv[1].split('=')[0] + '_eval.npy', out, allow_pickle=True)
