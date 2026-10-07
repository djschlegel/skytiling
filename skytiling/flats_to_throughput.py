#!/usr/bin/env python3
"""
Build a compact focal-plane throughput model from a set of LSST-pipeline flats.

Input: one afw ``flat`` file per (sub-filter, detector), e.g. the HSC MBQ1
quadrant-filter dome flats, which carry

* the normalised flat image (relative response),
* a MASK plane whose NO_DATA bit marks pixels not illuminated by that
  sub-filter (the obscuring cross, the vignetted field edge, dead amplifiers),
* the full camera geometry (see ``lsst_camera.py``).

Output: a single FITS file in field-angle coordinates (gnomonic offset from the
boresight, degrees, sampled at ``--res`` arcmin), with for each sub-filter

* ``GOOD_<f>``  fraction of unmasked (illuminated) pixels in each cell: 0 or 1
                away from edges, NaN where there is no detector at all,
* ``THRU_<f>``  relative throughput = normalised flat value x GOOD, 0 where
                masked, NaN where there is no detector.  The flat is normalised
                per sub-filter to 1 at its median over the inner field
                (r < ``--norm-radius``), so values are relative within a
                sub-filter only (dome-flat levels are not comparable between
                filters),

plus ``DETECTORS``/``FP2FA`` (camera geometry) and ``DETFLATS`` (per-file
statistics).  Field angle (x, y) is the LSST FIELD_ANGLE system of the camera.

Example::

    flats_to_throughput.py --flats-dir data/subaru/scratch/repo \\
        --out data/subaru/hsc_mbq1_throughput.fits --plot
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from skytiling.lsst_camera import LsstCamera

HSC_MBQ1_PATTERN = r'flat_mbq1_(?P<filter>\d+)_.*MBQ1_MBQ1_(?P<detector>\d_\d\d)_'


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--flats-dir', required=True, help='directory searched recursively for *.fits flats')
    p.add_argument('--pattern', default=HSC_MBQ1_PATTERN,
                   help='regex with named groups "filter" and "detector" applied to each file path')
    p.add_argument('--out', required=True, help='output FITS file')
    p.add_argument('--res', type=float, default=0.1, help='output cell size [arcmin]')
    p.add_argument('--extent', type=float, default=50.0, help='half-size of output grid [arcmin]')
    p.add_argument('--bin', type=int, default=8, help='input pixel binning before rasterising')
    p.add_argument('--norm-radius', type=float, default=20.0,
                   help='normalise each sub-filter flat to its median inside this field radius [arcmin]')
    p.add_argument('--mask-plane', default='NO_DATA', help='mask plane marking un-illuminated pixels')
    p.add_argument('--clip', type=float, default=2.0, help='clip flat values above this (hot pixels)')
    p.add_argument('--plot', action='store_true', help='write <out>.png with the maps')
    return p.parse_args(argv)


def main(argv=None):
    from astropy.io import fits
    opts = parse_arguments(argv)
    rx = re.compile(opts.pattern)
    files = []
    for f in sorted(glob.glob(os.path.join(opts.flats_dir, '**', '*.fits'), recursive=True)):
        m = rx.search(f)
        if m:
            files.append((m.group('filter'), m.group('detector'), f))
    if not files:
        sys.exit(f'no files matching pattern under {opts.flats_dir}')
    filters = sorted({f[0] for f in files})
    print(f'{len(files)} flats, sub-filters {filters}')

    cam = LsstCamera.from_fits(files[0][2])
    print(f'camera: {len(cam.detectors)} detectors, plate scale {cam.plate_scale_arcsec():.4f} arcsec/mm')

    res, half, B = opts.res, opts.extent, opts.bin
    N = int(round(2 * half / res))
    acc = {}      # (filter, kind) -> (sum, count) rasters
    stats = []
    for filt, det, f in files:
        with fits.open(f, memmap=True) as h:
            # afw MaskedImage layout: PRIMARY (no data), then IMAGE/MASK/VARIANCE extensions
            ext = {x.header.get('EXTTYPE'): x for x in h[1:] if x.header.get('EXTTYPE')}
            names = [x.name for x in h]
            img = np.array(ext['IMAGE'].data, dtype=np.float32)
            mask = np.array(ext['MASK'].data)
            bit = ext['MASK'].header[f'MP_{opts.mask_plane}']
            # per-file detector bbox (guide/focus CCDs differ from science CCDs)
            d = h[names.index('Detector')].data[0]
            nx = int(d['bbox_max_x'] - d['bbox_min_x'] + 1); ny = int(d['bbox_max_y'] - d['bbox_min_y'] + 1)
        cam.detectors[det]['nx'], cam.detectors[det]['ny'] = nx, ny
        bad = (mask & (1 << bit)) != 0
        good = ~bad
        img = np.minimum(np.nan_to_num(img, nan=0.0), opts.clip)
        ny_b, nx_b = ny // B, nx // B
        g = good[:ny_b * B, :nx_b * B].reshape(ny_b, B, nx_b, B).mean(axis=(1, 3))
        fsum = np.where(good, img, 0.0)[:ny_b * B, :nx_b * B].reshape(ny_b, B, nx_b, B).mean(axis=(1, 3))
        # mean flat over the *good* pixels of each bin, times the good fraction
        thr = np.where(g > 0, fsum / np.maximum(g, 1e-9), 0.0) * g
        py, px = np.mgrid[0:ny_b, 0:nx_b]
        px = px * B + (B - 1) / 2.0; py = py * B + (B - 1) / 2.0
        u, w = cam.pixel_to_field_angle(det, px, py)
        u *= 60.0; w *= 60.0
        iu = np.floor((u + half) / res).astype(int); iv = np.floor((w + half) / res).astype(int)
        ok = (iu >= 0) & (iu < N) & (iv >= 0) & (iv < N)
        for kind, arr in (('GOOD', g), ('THRU', thr)):
            s, c = acc.setdefault((filt, kind), (np.zeros((N, N)), np.zeros((N, N))))
            np.add.at(s, (iv[ok], iu[ok]), arr[ok]); np.add.at(c, (iv[ok], iu[ok]), 1.0)
        med = float(np.median(img[good])) if good.any() else np.nan
        stats.append((filt, det, float(good.mean()), med))
        print(f'  {filt} {det}: illuminated {good.mean():6.1%}  median flat {med:.4f}', flush=True)

    # rasters, normalisation
    yy, xx = np.mgrid[0:N, 0:N]
    r = np.hypot((xx + 0.5) * res - half, (yy + 0.5) * res - half)
    hdus = [fits.PrimaryHDU()]
    ph = hdus[0].header
    ph['NFILT'] = (len(filters), 'number of sub-filters')
    for i, f in enumerate(filters):
        ph[f'FILT{i}'] = f
    ph['RES'] = (res, 'cell size [arcmin]')
    ph['EXTENT'] = (half, 'grid half-size [arcmin]; x = (i+0.5)*RES-EXTENT')
    ph['NORMRAD'] = (opts.norm_radius, 'THRU normalised to median inside this radius [arcmin]')
    ph['MASKPLN'] = (opts.mask_plane, 'mask plane treated as un-illuminated')
    ph['FLATSDIR'] = os.path.basename(os.path.abspath(opts.flats_dir))
    ph.add_comment('Field-angle maps of relative throughput per sub-filter; see flats_to_throughput.py')
    ph.add_comment('x, y = LSST FIELD_ANGLE (gnomonic offset from boresight, arcmin), origin at centre')
    hdus += cam.to_hdus()
    stat_cols = [fits.Column(name='filter', format='8A', array=np.array([s[0] for s in stats])),
                 fits.Column(name='detector', format='8A', array=np.array([s[1] for s in stats])),
                 fits.Column(name='frac_illuminated', format='E', array=np.array([s[2] for s in stats])),
                 fits.Column(name='median_flat', format='E', array=np.array([s[3] for s in stats]))]
    hdus.append(fits.BinTableHDU.from_columns(stat_cols, name='DETFLATS'))
    maps = {}
    for filt in filters:
        s, c = acc[(filt, 'GOOD')]
        good = np.where(c > 0, s / np.maximum(c, 1), np.nan)
        s, c = acc[(filt, 'THRU')]
        thru = np.where(c > 0, s / np.maximum(c, 1), np.nan)
        inner = (r < opts.norm_radius) & (good > 0.99) & np.isfinite(thru)
        norm = float(np.median(thru[inner])) if inner.any() else 1.0
        thru = thru / norm
        maps[filt] = (good, thru)
        for kind, arr in (('GOOD', good), ('THRU', thru)):
            hdu = fits.CompImageHDU(arr.astype(np.float32), name=f'{kind}_{filt}', compression_type='RICE_1')
            hdu.header['FILTER'] = filt
            hdu.header['RES'] = (res, 'arcmin per cell')
            hdu.header['EXTENT'] = (half, 'arcmin')
            if kind == 'THRU':
                hdu.header['NORM'] = (norm, 'raw flat value mapped to 1.0')
            hdus.append(hdu)
        area = np.nansum(good) * res ** 2 / 3600.0
        print(f'MBQ {filt}: illuminated area {area:.4f} deg^2, THRU norm {norm:.4f}, '
              f'mean THRU over illuminated {np.nanmean(thru[good > 0.5]):.3f}')
    fits.HDUList(hdus).writeto(opts.out, overwrite=True)
    print(f'wrote {opts.out}')

    if opts.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, len(filters), figsize=(5.5 * len(filters), 11))
        for j, filt in enumerate(filters):
            good, thru = maps[filt]
            for i, (kind, arr, vmax, cmap) in enumerate((('illuminated fraction', good, 1, 'Greys_r'),
                                                         ('relative throughput', thru, 1.2, 'viridis'))):
                ax = axes[i, j]
                im = ax.imshow(arr, origin='lower', extent=[-half, half, -half, half], vmin=0, vmax=vmax, cmap=cmap)
                ax.set_title(f'{filt}: {kind}'); ax.set_xlabel('field angle x [arcmin]'); ax.set_ylabel('y [arcmin]')
                plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        png = os.path.splitext(opts.out)[0] + '.png'
        plt.savefig(png, dpi=80); print(f'wrote {png}')


if __name__ == '__main__':
    main()
