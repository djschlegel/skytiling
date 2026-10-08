#!/usr/bin/env python3
"""
Modify the filter-holder cross in a throughput-map file.

The MBQ1 flats mask everything within ``--old-half-width`` (5.2') of either
axis of the focal plane as un-illuminated.  If the opaque part of the cross is
really narrower (``--new-half-width``, e.g. 4.09' for an 8.18'-wide cross),
the band in between is detector area that does receive light, vignetted by
the filter holder.  This script takes a ``flats_to_throughput`` output file
and, for every map cell that has a detector, is currently masked, and lies in
that band (and outside the new cross), sets ``GOOD_<f>`` to 1 and ``THRU_<f>``
to a vignetting ramp: linear in the distance from the axis, from the map's own
value just outside the old mask edge (the flats already show the ramp
starting ~5' further out) down to ``--band-throughput`` at the opaque edge
(``--band-profile ramp``, the default; ``flat`` sets the whole band to
``--band-throughput``).  Cells inside the new cross stay masked; nothing else
changes (dead amplifiers, field edge, CCD gaps).

Example::

    adjust_cross -i data/subaru/hsc_mbq1_throughput.fits -o data/subaru2/hsc_mbq1_throughput.fits \\
        --old-half-width 5.2 --new-half-width 4.09 --band-throughput 0.8 --plot
"""
import argparse
import os
import sys

import numpy as np


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, help='flats_to_throughput FITS file')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--old-half-width', type=float, default=5.2, help='half-width of the cross masked in the flats [arcmin]')
    p.add_argument('--new-half-width', type=float, default=4.09, help='half-width of the opaque cross [arcmin]')
    p.add_argument('--band-throughput', type=float, default=0.8, help='relative throughput at the opaque edge of the cross')
    p.add_argument('--band-profile', choices=['ramp', 'flat'], default='ramp',
                   help='ramp: linear from the map value at the old mask edge to --band-throughput at the new edge; flat: constant')
    p.add_argument('--edge-width', type=float, default=0.5, help='ring outside the old edge used to read the local map value [arcmin]')
    p.add_argument('--plot', action='store_true', help='write <output>.png showing the change')
    return p.parse_args(argv)


def main(argv=None):
    from astropy.io import fits
    opts = parse_arguments(argv)
    h = fits.open(opts.input)
    ph = h[0].header
    res, half = float(ph['RES']), float(ph['EXTENT'])
    filters = [ph[f'FILT{i}'] for i in range(ph['NFILT'])]
    n = h[f'GOOD_{filters[0]}'].data.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    x = (xx + 0.5) * res - half; y = (yy + 0.5) * res - half
    r = np.hypot(x, y)
    dist = np.minimum(np.abs(x), np.abs(y))          # distance to the nearer axis
    band = (dist >= opts.new_half_width) & (dist < opts.old_half_width)
    ph['CROSSOLD'] = (opts.old_half_width, 'cross half-width masked in the flats [arcmin]')
    ph['CROSSNEW'] = (opts.new_half_width, 'opaque cross half-width adopted [arcmin]')
    ph['CROSSTHR'] = (opts.band_throughput, 'throughput in the band between them')
    ph.add_comment(f'adjust_cross: band {opts.new_half_width}-{opts.old_half_width} arcmin from the axes unmasked at throughput {opts.band_throughput}')
    ph['CROSSPRF'] = (opts.band_profile, 'band profile: ramp or flat')
    near_x = np.abs(x) <= np.abs(y)          # the x axis is the nearer one
    changed = {}
    for f in filters:
        good = h[f'GOOD_{f}'].data; thru = h[f'THRU_{f}'].data
        sel = band & np.isfinite(good) & (good < 0.5)
        good = good.copy(); thru = thru.copy()
        if opts.band_profile == 'flat':
            val = np.full(sel.sum(), opts.band_throughput)
        else:
            # local map value just outside the old mask edge, read along the row (near x axis)
            # or column (near y axis) through each band cell
            ring = (dist >= opts.old_half_width) & (dist < opts.old_half_width + opts.edge_width) & (good > 0.5) & np.isfinite(thru)
            edge = np.full((n, n), np.nan)
            # for cells nearest the x axis: same column (same ix), ring cells at |x| in the ring on the same side
            for axis_near, idx_same, coord in ((True, xx, x), (False, yy, y)):
                m = ring & (near_x == axis_near)
                for side in (-1, 1):
                    mm = m & (np.sign(coord) == side)
                    # average ring value per column (or row) index
                    sums = np.bincount(idx_same[mm], weights=thru[mm], minlength=n)
                    cnt = np.bincount(idx_same[mm], minlength=n)
                    prof = np.where(cnt > 0, sums / np.maximum(cnt, 1), np.nan)
                    tgt = sel & (near_x == axis_near) & (np.sign(coord) == side)
                    edge[tgt] = prof[idx_same[tgt]]
            e = edge[sel]
            fallback = np.nanmedian(thru[ring & (r < 35)]) if np.isfinite(thru[ring & (r < 35)]).any() else 1.0
            e = np.where(np.isfinite(e), e, fallback)
            frac = (opts.old_half_width - dist[sel]) / (opts.old_half_width - opts.new_half_width)   # 0 at old edge, 1 at new
            val = e + (opts.band_throughput - e) * frac
        good[sel] = 1.0; thru[sel] = val
        h[f'GOOD_{f}'].data = good.astype(np.float32); h[f'THRU_{f}'].data = thru.astype(np.float32)
        changed[f] = sel
        a_good = np.nansum(good) * (res / 60) ** 2; a_thru = np.nansum(thru) * (res / 60) ** 2
        print(f'{f}: {sel.sum()} cells unmasked ({sel.sum() * (res / 60) ** 2:.4f} deg^2); illuminated area now {a_good:.4f} deg^2, '
              f'throughput-weighted {a_thru:.4f} deg^2')
    h.writeto(opts.output, overwrite=True)
    print(f'wrote {opts.output}')
    if opts.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        nf = len(filters)
        zoom = 16.0
        fig, axes = plt.subplots(2, nf, figsize=(5.0 * nf, 10))
        for k, f in enumerate(filters):
            ax = axes[0, k]
            thru = h[f'THRU_{f}'].data
            im = ax.imshow(thru, origin='lower', extent=[-half, half, -half, half], vmin=0, vmax=1.1, cmap='viridis', interpolation='nearest')
            ax.contour(x, y, changed[f].astype(float), levels=[0.5], colors='r', linewidths=0.5)
            ax.set_xlim(-zoom, zoom); ax.set_ylim(-zoom, zoom); ax.set_aspect('equal')
            for v in (-opts.new_half_width, opts.new_half_width): ax.axhline(v, color='w', lw=0.4, ls=':'); ax.axvline(v, color='w', lw=0.4, ls=':')
            ax.set_title(f'MBQ1-{f}: central {2*zoom:.0f}\' at full resolution ({res}\'); red = ramp band', fontsize=9)
            ax.set_xlabel('field angle x [arcmin]'); ax.set_ylabel('y [arcmin]')
            plt.colorbar(im, ax=ax, fraction=0.046)
            ax = axes[1, k]
            good = h[f'GOOD_{f}'].data
            for axis_near, lab, col in ((True, 'distance from the y axis (|x|)', 'C0'), (False, 'distance from the x axis (|y|)', 'C1')):
                d = np.abs(x) if axis_near else np.abs(y)
                other = np.abs(y) if axis_near else np.abs(x)
                bins = np.arange(3.0, 16.0, 0.25); prof = []
                for lo in bins:
                    sl = (d >= lo) & (d < lo + 0.25) & (other > 8) & (r < 35) & (good > 0.99) & np.isfinite(thru)
                    prof.append(np.median(thru[sl]) if sl.sum() > 10 else np.nan)
                ax.plot(bins + 0.125, prof, '-', color=col, label=lab)
            ax.axvline(opts.new_half_width, color='k', lw=0.6, ls='--'); ax.axvline(opts.old_half_width, color='0.5', lw=0.6, ls=':')
            ax.axhline(opts.band_throughput, color='0.5', lw=0.5, ls=':')
            ax.set_ylim(0.6, 1.1); ax.set_xlim(3, 16); ax.grid(alpha=0.3)
            ax.set_xlabel('distance from the nearer axis [arcmin]'); ax.set_ylabel('median throughput')
            ax.set_title(f'{f}: profile toward the cross (dashed = opaque edge, dotted = old mask edge)', fontsize=9)
            if k == 0: ax.legend(fontsize=8, loc='lower right')
        plt.tight_layout()
        png = os.path.splitext(opts.output)[0] + '.png'
        plt.savefig(png, dpi=110); print(f'wrote {png}')


if __name__ == '__main__':
    main()
