#!/usr/bin/env python3
"""
Write the final tiling table (CSV + FITS) from an optimize_tiling_maps checkpoint, with
the cross-reference back to the starting layout: for a start built from pointing centres
x dither offsets, each row carries the pointing name (e.g. the HSC-SSP Wide pointing
``W_010001``), the dither index (1..n in the order of ``--init-offsets``), and the
starting (undrifted) position ``ra_init, dec_init`` of that dither, alongside the annealed
``ra, dec, rot``.  Positions are rounded to ``--ndec`` decimals and ``rot`` to 0.1 deg.

Example::

    write_tiling_table -i final_p2_tiles_0200.fits -o hsc_mbq1_tiling_stripes_final \\
        --offsets "12.2,19.5 -19.5,12.2 -12.2,-19.5 19.5,-12.2" --rotations 0,90,180,270 \\
        --comment "HSC MBQ1 tiling of the DESI Run 2 Niji stripes (data/subaru3)"
"""
import argparse

import numpy as np


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', required=True, help='tiles FITS checkpoint (tileid, ra, dec, rot, ra0, dec0[, pid, pname])')
    p.add_argument('-o', '--output', required=True, help='output base name (writes <base>.csv and <base>.fits)')
    p.add_argument('--offsets', help='the --init-offsets string used for the start, to label dithers [arcmin]')
    p.add_argument('--rotations', default='0,90,180,270', help='rotation set added to rot [deg]')
    p.add_argument('--ndec', type=int, default=3, help='decimals for ra, dec')
    p.add_argument('--comment', action='append', default=[], help='header comment line(s)')
    return p.parse_args(argv)


def main(argv=None):
    from astropy.io import fits
    opts = parse_arguments(argv)
    d = fits.getdata(opts.input)
    n = len(d)
    ra = np.round(d['ra'], opts.ndec); dec = np.round(d['dec'], opts.ndec); rot = np.round(d['rot'] % 360.0, 1)
    ra0 = np.round(d['ra0'], opts.ndec); dec0 = np.round(d['dec0'], opts.ndec)
    has_p = 'pid' in d.names
    pname = np.array(d['pname']) if 'pname' in d.names else (np.array([f'{p}' for p in d['pid']]) if has_p else np.array([''] * n))
    noff = len(opts.offsets.split()) if opts.offsets else 1
    dither = (np.arange(n) % noff) + 1 if has_p else np.zeros(n, int)
    cols = [fits.Column(name='tileid', format='K', array=d['tileid']),
            fits.Column(name='pointing', format=f'{max(1, max(len(s) for s in pname))}A', array=pname),
            fits.Column(name='dither', format='I', array=dither),
            fits.Column(name='ra', format='D', unit='deg', array=ra), fits.Column(name='dec', format='D', unit='deg', array=dec),
            fits.Column(name='rot', format='D', unit='deg', array=rot),
            fits.Column(name='ra_init', format='D', unit='deg', array=ra0), fits.Column(name='dec_init', format='D', unit='deg', array=dec0)]
    hdu = fits.BinTableHDU.from_columns(cols, name='TILES')
    hdu.header['ROTS'] = (opts.rotations, 'each tile is observed at rot + these [deg]')
    if opts.offsets:
        hdu.header['OFFSETS'] = (opts.offsets, 'dither offsets of the start [arcmin, +x = +RA]')
    for c in opts.comment:
        hdu.header.add_comment(c)
    for c in ('pointing, dither: pointing centre and dither index of the starting layout', 'that this tile descends from; ra_init, dec_init: that starting position;',
              'ra, dec, rot: the annealed position and rotation offset'):
        hdu.header.add_comment(c)
    hdu.writeto(f'{opts.output}.fits', overwrite=True)
    rots = ', '.join(f'rot+{r}' if float(r) else 'rot' for r in opts.rotations.split(','))
    with open(f'{opts.output}.csv', 'w') as f:
        for c in opts.comment:
            f.write(f'# {c}\n')
        f.write(f'# each row is one pointing centre, observed at rotation offsets {rots} deg\n')
        f.write('# pointing, dither: the pointing centre (HSC-SSP Wide name) and dither index of the starting layout this tile descends from\n')
        if opts.offsets:
            f.write(f'# dither offsets (arcmin, +x = +RA): {" | ".join(f"{i + 1}: ({o})" for i, o in enumerate(opts.offsets.split()))}\n')
        f.write('# ra_init, dec_init: that starting position; ra, dec, rot: the annealed position and rotation offset\n')
        f.write('tileid,pointing,dither,ra,dec,rot,ra_init,dec_init\n')
        for i in range(n):
            f.write(f'{d["tileid"][i]},{pname[i]},{dither[i]},{ra[i]:.{opts.ndec}f},{dec[i]:.{opts.ndec}f},{rot[i]:.1f},{ra0[i]:.{opts.ndec}f},{dec0[i]:.{opts.ndec}f}\n')
    print(f'wrote {opts.output}.csv and .fits ({n} rows)')


if __name__ == '__main__':
    main()
