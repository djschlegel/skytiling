#!/usr/bin/env python3
"""
Modify the filter-holder cross in a throughput-map file.

The MBQ1 flats mask everything within ``--old-half-width`` (5.2') of either
axis of the focal plane as un-illuminated.  If the opaque part of the cross is
really narrower (``--new-half-width``, e.g. 4.09' for an 8.18'-wide cross),
the band in between is detector area that does receive light, at a reduced
(vignetted) throughput ``--band-throughput``.  This script takes a
``flats_to_throughput`` output file and, for every map cell that has a
detector, is currently masked, and lies in that band (and outside the new
cross), sets ``GOOD_<f>`` to 1 and ``THRU_<f>`` to the band throughput.  Cells
inside the new cross stay masked; nothing else changes (dead amplifiers,
field edge, CCD gaps).

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
    p.add_argument('--band-throughput', type=float, default=0.8, help='relative throughput in the newly illuminated band')
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
    dist = np.minimum(np.abs(x), np.abs(y))          # distance to the nearer axis
    band = (dist >= opts.new_half_width) & (dist < opts.old_half_width)
    ph['CROSSOLD'] = (opts.old_half_width, 'cross half-width masked in the flats [arcmin]')
    ph['CROSSNEW'] = (opts.new_half_width, 'opaque cross half-width adopted [arcmin]')
    ph['CROSSTHR'] = (opts.band_throughput, 'throughput in the band between them')
    ph.add_comment(f'adjust_cross: band {opts.new_half_width}-{opts.old_half_width} arcmin from the axes unmasked at throughput {opts.band_throughput}')
    changed = {}
    for f in filters:
        good = h[f'GOOD_{f}'].data; thru = h[f'THRU_{f}'].data
        sel = band & np.isfinite(good) & (good < 0.5)
        good = good.copy(); thru = thru.copy()
        good[sel] = 1.0; thru[sel] = opts.band_throughput
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
        fig, axes = plt.subplots(1, len(filters), figsize=(5.2 * len(filters), 5.2))
        for ax, f in zip(axes, filters):
            thru = h[f'THRU_{f}'].data
            im = ax.imshow(thru, origin='lower', extent=[-half, half, -half, half], vmin=0, vmax=1.2, cmap='viridis')
            ax.contour(x, y, changed[f].astype(float), levels=[0.5], colors='r', linewidths=0.6)
            ax.set_title(f'MBQ1-{f}: throughput, cross {2 * opts.new_half_width:.2f}\' wide; red = band at {opts.band_throughput}')
            ax.set_xlabel('field angle x [arcmin]'); ax.set_ylabel('y [arcmin]')
            plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        png = os.path.splitext(opts.output)[0] + '.png'
        plt.savefig(png, dpi=80); print(f'wrote {png}')


if __name__ == '__main__':
    main()
