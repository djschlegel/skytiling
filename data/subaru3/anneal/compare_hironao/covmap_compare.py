import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/repo')
from skytiling.optimize_tiling_maps import *
from skytiling.optimize_tiling_maps import _tile_coverage
from astropy.io import fits
cam = ThroughputCamera('/home/claude/repo/data/subaru3/hsc_mbq1_throughput.fits')
def grid(tt, rot, ra0, dec0, diam, res=0.25):
    o = MapTilingOptimizer.__new__(MapTilingOptimizer); o.cam = cam; o.rot = rot
    o.tiles = radec_to_xyz(tt['ra'], tt['dec']); o.theta = np.radians(tt['rot'])
    o._tile_coverage = lambda t, pts, th=0.0: _tile_coverage(cam, rot, t, pts, th)
    return o.grid_coverage(ra0, dec0, diam, res)
filters = ['413', '439', '465', '490']
rot4 = [np.radians(r) for r in (0, 90, 180, 270)]
fig, axes = plt.subplots(2, 4, figsize=(17, 9.5))
for i, (f, rot, label) in enumerate((('/home/claude/hsc/stripes_eval/final_p2_tiles_0200.fits', rot4, 'skytiling stripes (rms/mean$^2$)'), ('/home/claude/hironao/wd/anneal_results/holes_rs4m5_l16_tiles_0300.fits', [0.0], 'Hironao rs4m5_l16 (var + 16 x holes)'))):
    tt = fits.getdata(f); X, Y, g = grid(tt, rot, 180.0, 1.5, 9.0)
    h = 4.5
    for k, fl in enumerate(filters):
        ax = axes[i, k]
        im = ax.imshow(g[:, :, k], origin='lower', extent=[h, -h, -h, h], vmin=0, vmax=6, cmap='viridis')
        for e in (-2.12 - 1.5, 0.68 - 1.5, 3.48 - 1.5): ax.axhline(e, color='r', lw=0.8)
        ax.set_title(f'MBQ1-{fl}: {label}', fontsize=10); ax.set_xlabel('dRA cos(Dec) [deg]'); ax.set_ylabel('dDec [deg] from Dec 1.5')
        plt.colorbar(im, ax=ax, fraction=0.046, label='depth' if k == 3 else None)
fig.suptitle('Throughput-weighted depth in a 9 x 9 deg field at RA 180 (red: the stripe edges Dec -2.12, 0.68, 3.48); Hironao\'s tiling also covers the Dec 4.33 row of Atsushi\'s pointings', fontsize=11)
plt.tight_layout(); plt.savefig('compare_covmap_ngc.png', dpi=80); print('done')
