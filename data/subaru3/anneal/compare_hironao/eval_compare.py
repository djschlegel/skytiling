"""Compare tilings on one set of independent randoms: depth (sum of THRU) and n_exp (sum of GOOD).

Regions: the DESI Run 2 Niji stripes (inside / interior = 0.75 deg inside) and Atsushi's 802-pointing
hex-cell footprint (inside / interior).  Usage: eval_compare.py out.npy name=tiles.fits ...
"""
import sys, numpy as np
sys.path.insert(0, '/home/claude/repo')
from skytiling.optimize_tiling_maps import *
from skytiling.optimize_tiling_maps import _tile_coverage
from astropy.io import fits
from scipy.spatial import cKDTree
import csv

class DualCamera(ThroughputCamera):
    def __init__(self, path):
        super().__init__(path, use_good=False)
        with fits.open(path) as h:
            good = np.stack([np.nan_to_num(h[f'GOOD_{f}'].data.astype(np.float32), nan=0.0) for f in self.filters])
        self.nthru = self.nf; self.maps = np.concatenate([self.maps, good]); self.nf = 2 * self.nthru

cam = DualCamera('/home/claude/repo/data/subaru3/hsc_mbq1_throughput.fits')
nt = cam.nthru
stripes = parse_rects('/home/claude/repo/data/subaru3/niji_stripes_2p8_hsc.csv')
rows = list(csv.DictReader(open('/home/claude/repo/data/subaru/hsc_niji_wide_pointings.csv')))
hexfp = Footprint(np.array([float(r['ra']) for r in rows]), np.array([float(r['dec']) for r in rows]), 0.75)
rng = np.random.default_rng(11)
N = 2_000_000
OUT = sys.argv[1]; args = sys.argv[2:]
# uniform randoms over the union of both footprints + 1.2 deg (rejection from two boxes)
def box_uniform(n, rng):
    out = []
    for lo, hi, dlo, dhi in ((120, 235, -5, 9), (320, 420, -9, 9)):
        m = n // 2
        rr = rng.uniform(lo, hi, m) % 360; z = rng.uniform(np.sin(np.radians(dlo)), np.sin(np.radians(dhi)), m)
        out.append(radec_to_xyz(rr, np.degrees(np.arcsin(z))))
    return np.vstack(out)
rand = box_uniform(int(N * 2.2), rng)
keep = stripes.contains(rand, 1.2) | hexfp.contains(rand, 1.2)
rand = rand[keep][:N]
ra, dec = xyz_to_radec(rand)
regions = dict(stripes_in=stripes.contains(rand), stripes_int=stripes.interior(rand, 0.75),
               hex_in=hexfp.contains(rand), hex_int=hexfp.interior(rand))
# per-rectangle masks
names = ['NGC-5', 'NGC-6', 'SGC-1', 'SGC-2', 'SGC-3', 'SGC-4']
for nm, (rlo, rhi, dlo, dhi) in zip(names, stripes.rects):
    regions['rect_' + nm] = (((ra - rlo) % 360) <= ((rhi - rlo) % 360)) & (dec >= dlo) & (dec <= dhi)
print('randoms', len(rand), {k: int(v.sum()) for k, v in regions.items()}, flush=True)
tree = cKDTree(rand); chord = 2 * np.sin(np.radians(cam.max_radius) / 2)
rot4 = [np.radians(r) for r in (0, 90, 180, 270)]; rot1 = [0.0]

def coverage(tt, rot):
    tiles = radec_to_xyz(tt['ra'], tt['dec']); theta = np.radians(tt['rot'])
    c = np.zeros((len(rand), cam.nf), np.float32)
    for i in range(len(tiles)):
        k = np.array(tree.query_ball_point(tiles[i], chord), dtype=np.int64)
        if len(k): c[k] += _tile_coverage(cam, rot, tiles[i], rand[k], theta[i])
    return c

def stats(c, sel):
    d = c[sel][:, :nt]; n = np.rint(c[sel][:, nt:])
    m, s = d.mean(0), d.std(0)
    return dict(mean=m, rms=s, cv=s / m, lt15=100 * (d < 1.5).mean(0), lt10=100 * (d < 1.0).mean(0),
                nle1=100 * (n <= 1).mean(0), nle1_union=100 * np.any(n <= 1, axis=1).mean(), nle0_union=100 * np.any(n <= 0, axis=1).mean(),
                lt15_union=100 * np.any(d < 1.5, axis=1).mean(), n_mean=n.mean(0))

out = {}
for a in args:
    name, path = a.split('=')
    tt = fits.getdata(path)
    rot = rot1 if (len(tt) > 8000 or 'ROTS' in fits.getheader(path, 1) and fits.getheader(path, 1)['ROTS'].strip() == '0') else rot4
    c = coverage(tt, rot)
    res = {k: stats(c, v) for k, v in regions.items()}
    tot = c[:, :nt].sum(0)
    res['spill_stripes'] = 100 * c[~regions['stripes_in']][:, :nt].sum(0) / tot
    res['spill_hex'] = 100 * c[~regions['hex_in']][:, :nt].sum(0) / tot
    tra, tdec = tt['ra'], tt['dec']; txyz = radec_to_xyz(tra, tdec)
    res['n_tiles'] = len(tt); res['n_exp'] = len(tt) * len(rot)
    res['n_exp_in_stripes'] = int(stripes.contains(txyz).sum()) * len(rot)
    res['n_exp_in_hex'] = int(hexfp.contains(txyz).sum()) * len(rot)
    # Dec profiles across the NGC band (RA 140-215) and the SGC band (RA 340-30)
    prof = {}
    for lab, rsel, dlo, dhi in (('NGC', (ra > 140) & (ra < 215), -3.5, 7.5), ('SGC', ((ra > 340) | (ra < 28)), -2.3, 7.5)):
        bins = np.arange(dlo, dhi + 0.001, 0.1); idx = np.digitize(dec, bins) - 1
        p = np.full((len(bins) - 1, cam.nf), np.nan)
        for b in range(len(bins) - 1):
            ss = rsel & (idx == b)
            if ss.sum() > 50: p[b] = c[ss].mean(0)
        prof[lab] = (0.5 * (bins[1:] + bins[:-1]), p)
    res['prof'] = prof
    res['hist'] = {k: np.histogram(c[regions[k]][:, 0], bins=np.arange(0, 8.01, 0.1))[0] for k in ('stripes_int', 'hex_int')}
    out[name] = res
    s = res['stripes_int']; h = res['hex_int']
    print(f"{name}: {res['n_exp']} exposures, {res['n_exp_in_stripes']} centred in the stripes; spill outside stripes {np.round(res['spill_stripes'],1)}%, outside hex {np.round(res['spill_hex'],1)}%")
    for lab, r in (('stripes interior', s), ('hex interior', h), ('stripes inside', res['stripes_in'])):
        print(f"   {lab:17s} mean {np.round(r['mean'],2)} cv {np.round(r['cv'],3)} sum {r['cv'].sum():.3f} | <1.5 {np.round(r['lt15'],1)}% union {r['lt15_union']:.1f}% | n<=1 {np.round(r['nle1'],1)}% union {r['nle1_union']:.2f}% | n=0 union {r['nle0_union']:.2f}%", flush=True)
    np.save(OUT, out, allow_pickle=True)
print('done')
