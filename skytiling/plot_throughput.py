#!/usr/bin/env python3
"""
Plot a throughput-map file (``flats_to_throughput`` / ``adjust_cross`` output),
one page per sub-filter: that filter's quadrant at full map resolution, the
throughput profiles toward the cross, and optionally the difference from a
second file (``--compare``), with the cells that are illuminated in only one of
the two outlined.

Example::

    plot_throughput -i data/subaru3/hsc_mbq1_throughput.fits --compare data/subaru2/hsc_mbq1_throughput.fits
"""
import argparse
import os

import numpy as np


def label(path):
    """Last two path components, e.g. subaru3/hsc_mbq1_throughput.fits."""
    return os.path.join(*os.path.normpath(os.path.abspath(path)).split(os.sep)[-2:])


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, help='throughput FITS file')
    p.add_argument('--compare', help='second throughput file; adds a difference panel (input minus this)')
    p.add_argument('-o', '--output', help='output base name (default: input without extension); writes <base>_<filter>.png')
    p.add_argument('--vmax', type=float, default=1.1)
    p.add_argument('--cross-half-width', type=float, default=4.09, help='opaque cross half-width to mark [arcmin]')
    return p.parse_args(argv)


def profiles(x, y, r, good, thru, step=0.25, lo=3.0, hi=20.0):
    out = {}
    for axis_near, key in ((True, 'x'), (False, 'y')):
        d = np.abs(x) if axis_near else np.abs(y)
        other = np.abs(y) if axis_near else np.abs(x)
        bins = np.arange(lo, hi, step); prof = []
        for b in bins:
            sl = (d >= b) & (d < b + step) & (other > 8) & (r < 35) & (good > 0.99) & np.isfinite(thru)
            prof.append(np.median(thru[sl]) if sl.sum() > 10 else np.nan)
        out[key] = (bins + step / 2, np.array(prof))
    return out


def main(argv=None):
    from astropy.io import fits
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    opts = parse_arguments(argv)
    h = fits.open(opts.input)
    c = fits.open(opts.compare) if opts.compare else None
    ph = h[0].header
    res, half = float(ph['RES']), float(ph['EXTENT'])
    filters = [ph[f'FILT{i}'] for i in range(ph['NFILT'])]
    n = h[f'GOOD_{filters[0]}'].data.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    x = (xx + 0.5) * res - half; y = (yy + 0.5) * res - half
    r = np.hypot(x, y)
    base = opts.output or os.path.splitext(opts.input)[0]
    cell = (res / 60) ** 2
    for f in filters:
        good = h[f'GOOD_{f}'].data; thru = h[f'THRU_{f}'].data
        ok = np.isfinite(good) & (good > 0.5)
        sx = 1 if np.nanmean(x[ok]) > 0 else -1; sy = 1 if np.nanmean(y[ok]) > 0 else -1
        ncol = 2 if c is not None else 1
        fig = plt.figure(figsize=(10 * ncol, 14))
        gs = fig.add_gridspec(2, ncol, height_ratios=[1.0, 0.42])
        ax = fig.add_subplot(gs[0, 0])
        im = ax.imshow(thru, origin='lower', extent=[-half, half, -half, half], vmin=0, vmax=opts.vmax, cmap='viridis', interpolation='nearest')
        ax.set_xlim(sorted([0.0, sx * half])); ax.set_ylim(sorted([0.0, sy * half])); ax.set_aspect('equal')
        ax.axhline(sy * opts.cross_half_width, color='w', lw=0.5, ls=':'); ax.axvline(sx * opts.cross_half_width, color='w', lw=0.5, ls=':')
        a_good = np.nansum(ok) * cell; a_thru = np.nansum(np.where(ok, thru, 0)) * cell
        ax.set_title(f'MBQ1-{f}: relative throughput at full map resolution ({res}\'); illuminated {a_good:.4f} deg$^2$, '
                     f'throughput-weighted {a_thru:.4f} deg$^2$\n{label(opts.input)}; dotted = opaque cross edge ({opts.cross_half_width}\')', fontsize=10)
        ax.set_xlabel('field angle x [arcmin]'); ax.set_ylabel('field angle y [arcmin]')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        ax = fig.add_subplot(gs[1, 0])
        pr = profiles(x, y, r, good, thru)
        ax.plot(*pr['x'], '-', color='C0', label='vs distance from the y axis (|x|)')
        ax.plot(*pr['y'], '-', color='C1', label='vs distance from the x axis (|y|)')
        if c is not None:
            cg = c[f'GOOD_{f}'].data; ct = c[f'THRU_{f}'].data
            pc = profiles(x, y, r, cg, ct)
            ax.plot(*pc['x'], ':', color='C0', label=f'|x|, {label(opts.compare)}')
            ax.plot(*pc['y'], ':', color='C1', label=f'|y|, {label(opts.compare)}')
        ax.axvline(opts.cross_half_width, color='k', lw=0.6, ls='--', label='opaque edge')
        ax.set_ylim(0.6, 1.1); ax.set_xlim(3, 20); ax.grid(alpha=0.3); ax.legend(fontsize=9, loc='lower right')
        ax.set_xlabel('distance from the nearer axis [arcmin]'); ax.set_ylabel('median throughput (r < 35\')')
        ax.set_title(f'MBQ1-{f}: throughput profile toward the cross', fontsize=10)
        if c is not None:
            cok = np.isfinite(cg) & (cg > 0.5)
            diff = np.where(ok & cok, thru - ct, np.nan)
            ax = fig.add_subplot(gs[0, 1])
            im = ax.imshow(diff, origin='lower', extent=[-half, half, -half, half], vmin=-0.2, vmax=0.2, cmap='RdBu_r', interpolation='nearest')
            ax.contour(x, y, (ok & ~cok).astype(float), levels=[0.5], colors='g', linewidths=0.7)
            ax.contour(x, y, (cok & ~ok).astype(float), levels=[0.5], colors='m', linewidths=0.7)
            ax.set_xlim(sorted([0.0, sx * half])); ax.set_ylim(sorted([0.0, sy * half])); ax.set_aspect('equal')
            new_a = np.sum(ok & ~cok) * cell; lost_a = np.sum(cok & ~ok) * cell
            ax.set_title(f'MBQ1-{f}: difference, {label(opts.input)} minus {label(opts.compare)}\n'
                         f'green outline: newly illuminated ({new_a:.4f} deg$^2$); magenta: no longer illuminated ({lost_a:.4f} deg$^2$)', fontsize=10)
            ax.set_xlabel('field angle x [arcmin]'); ax.set_ylabel('field angle y [arcmin]')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            ax = fig.add_subplot(gs[1, 1])
            rb = np.arange(0, 50, 0.5); med = []; medc = []
            for lo in rb:
                s = (r >= lo) & (r < lo + 0.5)
                med.append(np.nanmedian(thru[s & ok]) if (s & ok).sum() > 10 else np.nan)
                medc.append(np.nanmedian(ct[s & cok]) if (s & cok).sum() > 10 else np.nan)
            ax.plot(rb + 0.25, med, '-', label=label(opts.input))
            ax.plot(rb + 0.25, medc, ':', label=label(opts.compare))
            ax.set_ylim(0.6, 1.1); ax.set_xlim(0, 50); ax.grid(alpha=0.3); ax.legend(fontsize=9, loc='lower left')
            ax.set_xlabel('field radius [arcmin]'); ax.set_ylabel('median throughput of illuminated cells')
            ax.set_title(f'MBQ1-{f}: radial profile', fontsize=10)
        plt.tight_layout()
        png = f'{base}_{f}.png'
        plt.savefig(png, dpi=100); plt.close(fig); print(f'wrote {png}')


if __name__ == '__main__':
    main()
