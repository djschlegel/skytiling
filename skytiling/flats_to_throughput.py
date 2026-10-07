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
* ``THRU_<f>``  relative throughput = gain-corrected flat value x GOOD, 0 where
                masked, NaN where there is no detector.  The flats are in ADU,
                so each amplifier is multiplied by its gain (from the embedded
                Detector/Amplifier table) to get a response in electrons; the
                result is normalised per sub-filter to 1 at its median over the
                bulk of the field (r < ``--norm-radius``).  Values are relative
                within a sub-filter only (dome-flat levels are not comparable
                between filters).  Departures from 1 are the radial
                vignetting/illumination gradient, CCD-to-CCD QE differences and
                dead amplifiers,
                A flat from a uniform source also scales with the sky area per
                pixel, which changes across the field through the optical
                distortion (for HSC it falls to 0.90 of the on-axis value at
                r = 40' and 0.86 at 47'); a point source's flux does not, so
                that Jacobian (from the embedded distortion model) is divided
                out first.  The remaining smooth radial gradient (dome-screen
                illumination) is then divided out (``--flatten-radius``: a quadratic in field
                radius fitted to the azimuthal medians inside that radius and
                extrapolated beyond it), so the bulk of the field sits at ~1 and
                departures from 1 are the edge vignetting (steeper than the
                extrapolated trend), CCD-to-CCD QE differences, dead amplifiers
                and masked regions.  ``--flatten-radius 0`` keeps the gradient,
* ``RADIAL``    table of the azimuthally averaged throughput vs field radius
                per sub-filter (1' bins): the gain-corrected flat
                (``median_raw``), the sky area per pixel (``jacobian``), the
                fitted illumination trend (``model``) and the final throughput
                (``median_thru``),

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
    p.add_argument('--norm-radius', type=float, default=35.0,
                   help='normalise each sub-filter flat to its median inside this field radius [arcmin]')
    p.add_argument('--flatten-radius', type=float, default=38.0,
                   help='divide out a linear-in-radius fit to the azimuthal medians inside this '
                        'radius [arcmin]; 0 = do not flatten')
    p.add_argument('--flatten-rmin', type=float, default=12.0, help='inner radius of the trend fit [arcmin]')
    p.add_argument('--max-thru', type=float, default=1.0, help='cap THRU at this value (0 = no cap)')
    p.add_argument('--no-jacobian', action='store_true',
                   help='do not divide out the sky-area-per-pixel variation from the distortion model')
    p.add_argument('--no-gain', action='store_true',
                   help='do not multiply each amplifier by its gain (flats left in ADU)')
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
            det_hdus = [x for x in h if x.name == 'Detector']
            d = det_hdus[0].data[0]
            nx = int(d['bbox_max_x'] - d['bbox_min_x'] + 1); ny = int(d['bbox_max_y'] - d['bbox_min_y'] + 1)
            # second Detector HDU = amplifier table (bbox on the assembled image, gain)
            amps = [(int(a['bbox_min_x']), int(a['bbox_min_y']), int(a['bbox_extent_x']),
                     int(a['bbox_extent_y']), float(a['gain'])) for a in det_hdus[1].data]
        cam.detectors[det]['nx'], cam.detectors[det]['ny'] = nx, ny
        img = np.minimum(np.nan_to_num(img, nan=0.0), opts.clip)   # hot pixels, in the flat's own units
        if not opts.no_gain:
            for ax0, ay0, aw, ah, gain in amps:
                img[ay0:ay0 + ah, ax0:ax0 + aw] *= np.float32(gain)
        bad = (mask & (1 << bit)) != 0
        good = ~bad
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
        stats.append((filt, det, float(good.mean()), med, np.mean([a[4] for a in amps])))
        print(f'  {filt} {det}: illuminated {good.mean():6.1%}  median flat x gain {med:.4f}', flush=True)

    # rasters, normalisation
    yy, xx = np.mgrid[0:N, 0:N]
    xa = (xx + 0.5) * res - half; ya = (yy + 0.5) * res - half
    r = np.hypot(xa, ya)
    # Sky area per unit focal-plane area, relative to on-axis: the Jacobian of
    # the focal-plane -> field-angle map, evaluated on the output grid by
    # inverting the (nearly radial) map numerically.
    jac_map = np.ones((N, N))
    if not opts.no_jacobian:
        mm = np.linspace(0.0, 400.0, 4001)
        u_, v_ = cam.focal_plane_to_field_angle(mm, np.zeros_like(mm))
        r_of_mm = np.hypot(u_, v_) * 60.0
        keep = np.r_[True, np.diff(r_of_mm) > 0]
        mm_of_r = np.interp(r, r_of_mm[keep], mm[keep])
        with np.errstate(invalid='ignore', divide='ignore'):
            cx = np.where(r > 0, xa / r, 1.0); cy = np.where(r > 0, ya / r, 0.0)
        fx, fy = mm_of_r * cx, mm_of_r * cy
        d = 0.01
        u0, v0 = cam.focal_plane_to_field_angle(fx, fy)
        u1, v1 = cam.focal_plane_to_field_angle(fx + d, fy)
        u2, v2 = cam.focal_plane_to_field_angle(fx, fy + d)
        jac_map = np.abs((u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0)) / d ** 2
        j0 = jac_map[N // 2, N // 2]
        jac_map = jac_map / j0
        row = r[N // 2, N // 2:]; jrow = jac_map[N // 2, N // 2:]
        print('sky area per pixel relative to on-axis: ' +
              ', '.join(f'r={rq}\' {np.interp(rq, row, jrow):.3f}' for rq in (20, 30, 40, 47)))
    hdus = [fits.PrimaryHDU()]
    ph = hdus[0].header
    ph['NFILT'] = (len(filters), 'number of sub-filters')
    for i, f in enumerate(filters):
        ph[f'FILT{i}'] = f
    ph['RES'] = (res, 'cell size [arcmin]')
    ph['EXTENT'] = (half, 'grid half-size [arcmin]; x = (i+0.5)*RES-EXTENT')
    ph['NORMRAD'] = (opts.norm_radius, 'THRU normalised to median inside this radius [arcmin]')
    ph['GAINCORR'] = (not opts.no_gain, 'flats multiplied by per-amplifier gain')
    ph['JACCORR'] = (not opts.no_jacobian, 'sky area per pixel (distortion Jacobian) divided out')
    ph['FLATRAD'] = (opts.flatten_radius, 'radial trend fitted inside this radius and divided out (0=no)')
    ph['FLATRMIN'] = (opts.flatten_rmin, 'inner radius of the trend fit [arcmin]')
    ph['MAXTHRU'] = (opts.max_thru, 'THRU capped at this value (0=no cap)')
    ph['MASKPLN'] = (opts.mask_plane, 'mask plane treated as un-illuminated')
    ph['FLATSDIR'] = os.path.basename(os.path.abspath(opts.flats_dir))
    ph.add_comment('Field-angle maps of relative throughput per sub-filter; see flats_to_throughput.py')
    ph.add_comment('x, y = LSST FIELD_ANGLE (gnomonic offset from boresight, arcmin), origin at centre')
    hdus += cam.to_hdus()
    stat_cols = [fits.Column(name='filter', format='8A', array=np.array([s[0] for s in stats])),
                 fits.Column(name='detector', format='8A', array=np.array([s[1] for s in stats])),
                 fits.Column(name='frac_illuminated', format='E', array=np.array([s[2] for s in stats])),
                 fits.Column(name='median_flat', format='E', array=np.array([s[3] for s in stats])),
                 fits.Column(name='mean_gain', format='E', array=np.array([s[4] for s in stats]))]
    hdus.append(fits.BinTableHDU.from_columns(stat_cols, name='DETFLATS'))
    maps = {}
    radial = []
    rbins = np.arange(0.0, half, 1.0)
    for filt in filters:
        s, c = acc[(filt, 'GOOD')]
        good = np.where(c > 0, s / np.maximum(c, 1), np.nan)
        s, c = acc[(filt, 'THRU')]
        thru = np.where(c > 0, s / np.maximum(c, 1), np.nan) / jac_map
        inner = (r < opts.norm_radius) & (good > 0.99) & np.isfinite(thru)
        norm = float(np.median(thru[inner])) if inner.any() else 1.0
        thru = thru / norm
        # azimuthal medians of the gain- and Jacobian-corrected, normalised flat
        raw_prof = []
        for rb in rbins:
            sel = (r >= rb) & (r < rb + 1) & (good > 0.99) & np.isfinite(thru)
            jsel = (r >= rb) & (r < rb + 1)
            raw_prof.append((rb + 0.5, float(np.median(thru[sel])) if sel.sum() > 20 else np.nan, int(sel.sum()),
                             float(np.median(jac_map[jsel]))))
        rp = np.array([(a, b) for a, b, c, j in raw_prof if np.isfinite(b)])
        if opts.flatten_radius > 0:
            fit = rp[(rp[:, 0] < opts.flatten_radius) & (rp[:, 0] > opts.flatten_rmin)]
            coef = np.polyfit(fit[:, 0], fit[:, 1], 1)
            model = np.polyval(coef, r)
            thru = thru / model
            if opts.max_thru > 0:
                thru = np.minimum(thru, opts.max_thru)
            model_prof = np.polyval(coef, np.array([a for a, b, c, j in raw_prof]))
        else:
            coef = np.array([0.0, 1.0]); model_prof = np.ones(len(raw_prof))
        maps[filt] = (good, thru, coef)
        for (rc, med_raw, n, jm), mdl in zip(raw_prof, model_prof):
            sel = (r >= rc - 0.5) & (r < rc + 0.5) & (good > 0.99) & np.isfinite(thru)
            radial.append((filt, rc, med_raw, jm, float(mdl), float(np.median(thru[sel])) if n > 20 else np.nan, n))
        for kind, arr in (('GOOD', good), ('THRU', thru)):
            hdu = fits.CompImageHDU(arr.astype(np.float32), name=f'{kind}_{filt}', compression_type='RICE_1')
            hdu.header['FILTER'] = filt
            hdu.header['RES'] = (res, 'arcmin per cell')
            hdu.header['EXTENT'] = (half, 'arcmin')
            if kind == 'THRU':
                hdu.header['NORM'] = (norm, 'gain-corrected flat value mapped to 1.0 before flattening')
                for k, c in enumerate(coef[::-1]):
                    hdu.header[f'RADC{k}'] = (c, f'radial trend coef of r^{k} [r in arcmin], divided out')
            hdus.append(hdu)
        area = np.nansum(good) * res ** 2 / 3600.0
        print(f'MBQ {filt}: illuminated area {area:.4f} deg^2, THRU norm {norm:.4f}, '
              f'mean THRU over illuminated {np.nanmean(thru[good > 0.5]):.3f}')
    rad_cols = [fits.Column(name='filter', format='8A', array=np.array([x[0] for x in radial])),
                fits.Column(name='radius', format='E', unit='arcmin', array=np.array([x[1] for x in radial])),
                fits.Column(name='median_raw', format='E', array=np.array([x[2] for x in radial])),
                fits.Column(name='jacobian', format='E', array=np.array([x[3] for x in radial])),
                fits.Column(name='model', format='E', array=np.array([x[4] for x in radial])),
                fits.Column(name='median_thru', format='E', array=np.array([x[5] for x in radial])),
                fits.Column(name='ncells', format='J', array=np.array([x[6] for x in radial]))]
    rad_hdu = hdus[-1] if False else None
    rad_hdu = fits.BinTableHDU.from_columns(rad_cols, name='RADIAL')
    rad_hdu.header['COMMENT'] = 'median_raw: gain-corrected flat / jacobian, normalised at r < NORMRAD'
    rad_hdu.header['COMMENT'] = 'jacobian: sky area per pixel relative to on-axis; model: fitted illumination trend'
    hdus.append(rad_hdu)
    fits.HDUList(hdus).writeto(opts.out, overwrite=True)
    print(f'wrote {opts.out}')

    if opts.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        nf = len(filters)
        fig = plt.figure(figsize=(5.5 * nf, 16.5))
        gs = fig.add_gridspec(3, nf, height_ratios=[1, 1, 0.8])
        for j, filt in enumerate(filters):
            good, thru, _ = maps[filt]
            for i, (kind, arr) in enumerate((('illuminated fraction', good), ('relative throughput', thru))):
                ax = fig.add_subplot(gs[i, j])
                im = ax.imshow(arr, origin='lower', extent=[-half, half, -half, half], vmin=0, vmax=1.2, cmap='viridis')
                ax.set_title(f'MBQ1-{filt}: {kind}'); ax.set_xlabel('field angle x [arcmin]'); ax.set_ylabel('y [arcmin]')
                plt.colorbar(im, ax=ax, fraction=0.046)
        ax = fig.add_subplot(gs[2, :])
        for k, filt in enumerate(filters):
            sel = [x for x in radial if x[0] == filt]
            rr = np.array([x[1] for x in sel]); raw = np.array([x[2] for x in sel]); jj = np.array([x[3] for x in sel])
            mdl = np.array([x[4] for x in sel]); fl = np.array([x[5] for x in sel])
            col = f'C{k}'
            if k == 0:
                ax.plot(rr, jj, '-', color='0.4', lw=1.5, label='sky area per pixel (distortion Jacobian)')
            ax.plot(rr, raw, ':', color=col, lw=1, label='gain-corrected flat / Jacobian' if k == 0 else None)
            ax.plot(rr, mdl, '--', color=col, lw=0.8, label='fitted illumination trend' if k == 0 else None)
            ax.plot(rr, fl, '-o', color=col, ms=3, label=f'MBQ1-{filt}: throughput')
        ax.axhline(1.0, color='k', lw=0.5); ax.axvline(opts.flatten_radius, color='k', lw=0.5, ls='--')
        ax.set_xlabel('field radius [arcmin]'); ax.set_ylabel('azimuthal median'); ax.set_ylim(0.6, 1.2)
        ax.set_title(f'radial profiles: flat / Jacobian normalised at r < {opts.norm_radius:g}\' (dotted), '
                     f'linear trend fitted at {opts.flatten_rmin:g}-{opts.flatten_radius:g}\' (dashed), '
                     f'throughput = min(flat / trend, {opts.max_thru:g}) (solid)')
        ax.legend(ncol=3, loc='lower left'); ax.grid(alpha=0.3)
        plt.tight_layout()
        png = os.path.splitext(opts.out)[0] + '.png'
        plt.savefig(png, dpi=80); print(f'wrote {png}')


if __name__ == '__main__':
    main()
