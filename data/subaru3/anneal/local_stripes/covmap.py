import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/repo')
from skytiling.optimize_tiling_maps import *
from skytiling.optimize_tiling_maps import _tile_coverage
from astropy.io import fits
cam = ThroughputCamera('/home/claude/repo/data/subaru3/hsc_mbq1_throughput.fits')
rot = [np.radians(r) for r in (0, 90, 180, 270)]
def grid(tt, ra0, dec0, diam, res=0.25):
    o = MapTilingOptimizer.__new__(MapTilingOptimizer); o.cam = cam; o.rot = rot
    o.tiles = radec_to_xyz(tt['ra'], tt['dec']); o.theta = np.radians(tt['rot'])
    o._tile_coverage = lambda t, pts, th=0.0: _tile_coverage(cam, rot, t, pts, th)
    return o.grid_coverage(ra0, dec0, diam, res)
filters = ['413', '439', '465', '490']
fig, axes = plt.subplots(2, 4, figsize=(17, 9))
for i, (f, label) in enumerate((('ssp_p1_tiles_0000.fits', 'start: SSP pointings + 4 dithers'), ('final_p2_tiles_0200.fits', 'annealed, rms/mean$^2$ objective'))):
    tt = fits.getdata(f); X, Y, g = grid(tt, 180.0, 0.68, 7.0)
    h = 3.5
    for k, fl in enumerate(filters):
        ax = axes[i, k]
        im = ax.imshow(g[:, :, k], origin='lower', extent=[h, -h, -h, h], vmin=0, vmax=6, cmap='viridis')
        for e in (-2.12 - 0.68, 3.48 - 0.68): ax.axhline(e, color='r', lw=0.8)
        ax.set_title(f'MBQ1-{fl}: {label}', fontsize=10); ax.set_xlabel('dRA cos(Dec) [deg]'); ax.set_ylabel('dDec [deg] from Dec 0.68')
        plt.colorbar(im, ax=ax, fraction=0.046, label='coverage' if k == 3 else None)
fig.suptitle('Throughput-weighted coverage in a 7 x 7 deg field at RA 180 across the NGC band (red: the band edges, Dec -2.12 and 3.48); same number of exposures (2720 pointings x 4 rotations)', fontsize=11)
plt.tight_layout(); plt.savefig('hsc_stripes_coverage_compare.png', dpi=90); print('done')
