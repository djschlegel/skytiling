#!/usr/bin/env python3
"""
Evaluate the coverage of a regular-grid + dither tiling with per-sub-filter
throughput maps (as written by ``flats_to_throughput``).

Geometry follows hsc-niji-survey/wide-observing-strategy/wide_dithering: a
hex grid of pointings (``--col-spacing`` x ``--row-spacing`` arcmin), every
pointing observed at the same set of dither offsets (arcmin, +x = +RA,
+y = +Dec), each dither at every instrument rotation in ``--rotations``.
Coverage is accumulated on one periodic cell of the grid (centred on a
pointing) from every pointing within ``--reach`` of the cell, which is exact
for dither radii up to reach - 55'.

For each sub-filter two quantities are accumulated per sky cell:

* n_obs  = number of exposures whose illuminated footprint (``GOOD_<f>``)
           covers the cell (0/1 coverage, as in Atsushi's evaluation);
* depth  = sum over exposures of the relative throughput (``THRU_<f>``),
           i.e. the effective number of full-throughput exposures.

Reported per sub-filter: mean, RMS, RMS/mean, fraction with n_obs = 0 and
<= 1, and the same for depth with thresholds 0.5 and 1.5.

Examples::

    eval_dither_tiling -t data/subaru/hsc_mbq1_throughput.fits --pattern 4plus2
    eval_dither_tiling -t ... --offsets "0,0 8,29.9 -8,-29.9"
"""
import argparse
import os
import sys

import numpy as np

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Atsushi's recommended patterns (wide_dithering/README.md), offsets in arcmin (dRA, dDec)
PATTERNS = {
    'none': [(0.0, 0.0)],
    '4plus2_phase1': [(12.2, 19.5), (-19.5, 12.2), (-12.2, -19.5), (19.5, -12.2)],
    '4plus2': [(12.2, 19.5), (-19.5, 12.2), (-12.2, -19.5), (19.5, -12.2), (30.3, 17.5), (-30.3, -17.5)],
    '4plus4': [(12.2, 19.5), (-19.5, 12.2), (-12.2, -19.5), (19.5, -12.2)],   # phase 2 square not quoted numerically
    '3plus3_phase1': [(0.0, 0.0), (8.0, 29.9), (-8.0, -29.9)],
    '3plus3': [(0.0, 0.0), (8.0, 29.9), (-8.0, -29.9), (29.1, 7.3), (-29.1, -7.3), (20.8, -21.6)],
}


def hex_pointings_near_cell(reach, col, row):
    """Hex-grid pointing centres (arcmin) within `reach` of the central cell."""
    n = 2 * int(np.ceil(reach / min(col, row))) + 3
    half = n // 2
    pts = []
    for j in range(-half, half + 1):
        xoff = col / 2.0 if j % 2 else 0.0
        for i in range(-half, half + 1):
            pts.append((i * col + xoff, j * row))
    p = np.array(pts)
    dx = np.maximum(np.abs(p[:, 0]) - col / 2, 0.0)
    dy = np.maximum(np.abs(p[:, 1]) - row / 2, 0.0)
    return p[np.hypot(dx, dy) <= reach]


class ThroughputMaps:
    """Nearest-cell lookup into the field-angle maps of a flats_to_throughput file."""

    def __init__(self, path):
        from astropy.io import fits
        with fits.open(path) as h:
            ph = h[0].header
            self.filters = [ph[f'FILT{i}'] for i in range(ph['NFILT'])]
            self.res = float(ph['RES']); self.half = float(ph['EXTENT'])
            self.good = {f: np.nan_to_num(h[f'GOOD_{f}'].data.astype(np.float32), nan=0.0) for f in self.filters}
            self.thru = {f: np.nan_to_num(h[f'THRU_{f}'].data.astype(np.float32), nan=0.0) for f in self.filters}
        self.n = self.good[self.filters[0]].shape[0]

    def lookup(self, arr, lx, ly):
        ix = np.floor((lx + self.half) / self.res).astype(np.int64)
        iy = np.floor((ly + self.half) / self.res).astype(np.int64)
        ok = (ix >= 0) & (ix < self.n) & (iy >= 0) & (iy < self.n)
        out = np.zeros(lx.shape, dtype=np.float32)
        out[ok] = arr[iy[ok], ix[ok]]
        return out


def coverage_maps(maps, centers, rotations_deg, res, col, row):
    """Return x, y grids and dicts n_obs[f], depth[f] over the periodic cell."""
    xg = np.arange(-col / 2, col / 2, res) + res / 2
    yg = np.arange(-row / 2, row / 2, res) + res / 2
    X, Y = np.meshgrid(xg, yg)
    nobs = {f: np.zeros(X.shape, np.float32) for f in maps.filters}
    depth = {f: np.zeros(X.shape, np.float32) for f in maps.filters}
    for phi_deg in rotations_deg:
        c, s = np.cos(np.radians(phi_deg)), np.sin(np.radians(phi_deg))
        for cx, cy in centers:
            dx, dy = X - cx, Y - cy
            if np.min(np.hypot(dx, dy)) > maps.half * 1.42:
                continue
            lx, ly = dx * c + dy * s, -dx * s + dy * c
            for f in maps.filters:
                nobs[f] += maps.lookup(maps.good[f], lx, ly) > 0.5
                depth[f] += maps.lookup(maps.thru[f], lx, ly)
    return xg, yg, nobs, depth


def summarize(nobs, depth, filters):
    rows = []
    for f in filters:
        n, d = nobs[f], depth[f]
        rows.append(dict(filter=f,
                         n_mean=n.mean(), n_rms=n.std(), n_cv=n.std() / n.mean(),
                         n0=(n == 0).mean(), n_le1=(n <= 1).mean(),
                         d_mean=d.mean(), d_rms=d.std(), d_cv=d.std() / d.mean(),
                         d_lt05=(d < 0.5).mean(), d_lt15=(d < 1.5).mean()))
    return rows


def print_table(rows, label):
    print(f'\n{label}')
    print(f'{"filter":>7} | {"n_obs mean":>10} {"rms":>6} {"rms/mean":>8} {"n=0":>7} {"n<=1":>7} | '
          f'{"depth mean":>10} {"rms":>6} {"rms/mean":>8} {"d<0.5":>7} {"d<1.5":>7}')
    for r in rows:
        print(f'{r["filter"]:>7} | {r["n_mean"]:10.3f} {r["n_rms"]:6.3f} {r["n_cv"]:8.3f} '
              f'{100 * r["n0"]:6.2f}% {100 * r["n_le1"]:6.2f}% | '
              f'{r["d_mean"]:10.3f} {r["d_rms"]:6.3f} {r["d_cv"]:8.3f} '
              f'{100 * r["d_lt05"]:6.2f}% {100 * r["d_lt15"]:6.2f}%')
    m = lambda k: np.mean([r[k] for r in rows])
    print(f'{"mean":>7} | {m("n_mean"):10.3f} {m("n_rms"):6.3f} {m("n_cv"):8.3f} '
          f'{100 * m("n0"):6.2f}% {100 * m("n_le1"):6.2f}% | '
          f'{m("d_mean"):10.3f} {m("d_rms"):6.3f} {m("d_cv"):8.3f} '
          f'{100 * m("d_lt05"):6.2f}% {100 * m("d_lt15"):6.2f}%')


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-t', '--throughput', required=True, help='flats_to_throughput FITS file')
    p.add_argument('--pattern', default='4plus2', help='named dither pattern: ' + ', '.join(PATTERNS))
    p.add_argument('--offsets', help='explicit dither offsets "dx,dy dx,dy ..." [arcmin], overrides --pattern')
    p.add_argument('--col-spacing', type=float, default=77.5, help='hex grid column spacing [arcmin]')
    p.add_argument('--row-spacing', type=float, default=67.5, help='hex grid row spacing [arcmin]')
    p.add_argument('--rotations', default='0,90,180,270', help='instrument rotations [deg]')
    p.add_argument('--res', type=float, default=0.1, help='evaluation cell size [arcmin]')
    p.add_argument('--reach', type=float, default=120.0, help='include pointings within this distance of the cell')
    p.add_argument('--plot', help='write coverage maps to this PNG')
    return p.parse_args(argv)


def main(argv=None):
    opts = parse_arguments(argv)
    maps = ThroughputMaps(opts.throughput)
    if opts.offsets:
        offsets = np.array([[float(v) for v in o.split(',')] for o in opts.offsets.split()])
        label = 'custom'
    else:
        offsets = np.array(PATTERNS[opts.pattern]); label = opts.pattern
    rot = [float(v) for v in opts.rotations.split(',')]
    pointings = hex_pointings_near_cell(opts.reach, opts.col_spacing, opts.row_spacing)
    centers = (pointings[:, None, :] + offsets[None, :, :]).reshape(-1, 2)
    cell_area = opts.col_spacing * opts.row_spacing / 3600.0
    print(f'{len(pointings)} pointings near cell, {len(offsets)} dithers x {len(rot)} rotations '
          f'= {len(offsets) * len(rot)} exposures per pointing; cell {opts.col_spacing}\' x {opts.row_spacing}\' '
          f'= {cell_area:.3f} deg^2')
    for f in maps.filters:
        a_good = maps.good[f].sum() * maps.res ** 2 / 3600.0
        a_thru = maps.thru[f].sum() * maps.res ** 2 / 3600.0
        print(f'  {f}: illuminated {a_good:.4f} deg^2, throughput-weighted {a_thru:.4f} deg^2 -> expected mean '
              f'n_obs {len(offsets) * len(rot) * a_good / cell_area:.3f}, depth {len(offsets) * len(rot) * a_thru / cell_area:.3f}')
    xg, yg, nobs, depth = coverage_maps(maps, centers, rot, opts.res, opts.col_spacing, opts.row_spacing)
    rows = summarize(nobs, depth, maps.filters)
    print_table(rows, f'pattern {label}: {len(offsets)} dithers, rotations {rot}, {opts.res}\' cells')
    union_le1 = np.mean(np.any([nobs[f] <= 1 for f in maps.filters], axis=0))
    union_0 = np.mean(np.any([nobs[f] == 0 for f in maps.filters], axis=0))
    print(f'union over filters: n<=1 in >=1 filter {100 * union_le1:.2f}%, n=0 in >=1 filter {100 * union_0:.3f}%')

    if opts.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        nf = len(maps.filters)
        fig, axes = plt.subplots(2, nf, figsize=(4.5 * nf, 8.5))
        ext = [xg[0], xg[-1], yg[0], yg[-1]]
        vmax = max(np.percentile(depth[f], 99.5) for f in maps.filters) * 1.05
        for j, f in enumerate(maps.filters):
            for i, (kind, arr) in enumerate((('n_obs', nobs[f]), ('depth (throughput-weighted)', depth[f]))):
                ax = axes[i, j]
                im = ax.imshow(arr, origin='lower', extent=ext, vmin=0, vmax=vmax, cmap='viridis')
                ax.set_title(f'MBQ1-{f}: {kind}  mean {arr.mean():.2f}, rms {arr.std():.2f}', fontsize=9)
                ax.set_xlabel('dRA [arcmin]'); ax.set_ylabel('dDec [arcmin]')
                plt.colorbar(im, ax=ax, fraction=0.046)
        fig.suptitle(f'{label}: {len(offsets)} dithers x {len(rot)} rotations per pointing, hex grid '
                     f'{opts.col_spacing}\' x {opts.row_spacing}\'')
        plt.tight_layout(); plt.savefig(opts.plot, dpi=80); print(f'wrote {opts.plot}')


if __name__ == '__main__':
    main()
