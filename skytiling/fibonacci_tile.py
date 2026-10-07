#!/usr/bin/env python
"""
Generates an all-sky grid of points using a Fibonacci Sphere algorithm.
This provides a highly uniform, equal-area distribution of points across the
entire celestial sphere for any arbitrary point density.

Can be used as a command-line tool::

    fibonacci_tile.py -a 0.25 [-o basename]

or imported::

    from skytiling.fibonacci_tile import fibonacci_tiling
    result = fibonacci_tiling(0.25)
"""
import sys
import numpy as np
import argparse

__version__ = "1.1-fibonacci"

# Total surface area of a sphere in square degrees: 4 pi (180/pi)^2 ~ 41252.9612
SPHERE_AREA_DEG2 = 4 * np.pi * (180 / np.pi)**2


def fibonacci_points(ntiles):
    """Return (tileid, ra, dec) arrays for an ntiles-point Fibonacci lattice.

    RA in [0, 360), Dec in [-90, 90], both in degrees.  tileid is 1-based.
    """
    # Spread points uniformly along a golden spiral wrapped around the sphere
    indices = np.arange(0, ntiles, dtype=float) + 0.5

    # Co-latitude angle phi (from 0 at the north pole to pi at the south pole)
    phi = np.arccos(1 - 2 * indices / ntiles)

    # Longitude angle theta using the Golden Ratio increment
    golden_ratio = (1 + 5**0.5) / 2
    theta = 2 * np.pi * golden_ratio * indices

    # Convert spherical coordinates to Celestial coordinates (RA and Dec)
    ra = np.degrees(theta) % 360.0
    dec = np.degrees(np.pi / 2.0 - phi)

    tileid = np.arange(ntiles) + 1
    return tileid, ra, dec


def fibonacci_tiling(tilearea, output=None, write_fits=True, write_ascii=True,
                     verbose=True):
    """Generate an all-sky Fibonacci tiling with ``tilearea`` deg^2 per tile.

    Parameters
    ----------
    tilearea : float
        Desired average area per tile in square degrees (e.g. 0.25).
    output : str, optional
        Output basename.  Defaults to ``allsky-fib-{ntiles}-{area:.3f}``.
    write_fits, write_ascii : bool
        Which output formats to write (``.fits`` requires astropy).
    verbose : bool
        Print progress to stdout.

    Returns
    -------
    dict with keys ``ntiles``, ``tilearea`` (actual), ``output`` (basename),
    ``tileid``, ``ra``, ``dec`` and ``files`` (list of files written).
    """
    if tilearea <= 0:
        raise ValueError("Tile area must be a positive number.")

    # Calculate the exact number of points (N) needed to satisfy the requested area
    ntiles = int(round(SPHERE_AREA_DEG2 / tilearea))

    # Recalculate actual average area per point due to integer rounding of N
    actual_tilearea = SPHERE_AREA_DEG2 / ntiles

    if verbose:
        print(f"Requested target area: {tilearea} sq deg/point")
        print(f"Generating {ntiles} points (Actual spacing: {actual_tilearea:.5f} sq deg/point)...")

    if output is None:
        output = f'allsky-fib-{ntiles}-{actual_tilearea:.3f}'

    tileid, ra, dec = fibonacci_points(ntiles)
    files = []

    # ----- Output formats
    if write_fits:
        try:
            from astropy.io import fits

            col1 = fits.Column(name='tileid', format='K', array=tileid)
            col2 = fits.Column(name='ra', format='D', array=ra)
            col3 = fits.Column(name='dec', format='D', array=dec)
            coldefs = fits.ColDefs([col1, col2, col3])

            hdu = fits.BinTableHDU.from_columns(coldefs, name='ALLSKY')

            hdu.header['NTILES'] = (ntiles, 'Number of tiles')
            hdu.header['TILEAREA'] = (actual_tilearea, 'Typical tile area [deg^2]')
            hdu.header.add_comment('Generated using analytical Fibonacci Lattice distribution')
            hdu.header.add_comment(f'Written by fibonacci_tile.py version {__version__}')
            hdu.header.add_comment('(ra,dec) in degrees')

            hdu.writeto(output + '.fits', overwrite=True)
            files.append(output + '.fits')
            if verbose:
                print(f"Wrote {output}.fits")
        except ImportError:
            print("NOTE: astropy package not installed. Skipping .fits creation.", file=sys.stderr)

    if write_ascii:
        with open(output + '.dat', 'w') as fx:
            print(f"# {ntiles} tiles, typical area {actual_tilearea:.5f} deg^2 per tile", file=fx)
            print("# Generated using analytical Fibonacci Lattice distribution", file=fx)
            print(f'# Written by fibonacci_tile.py version {__version__}', file=fx)
            print("# (ra,dec) in degrees", file=fx)
            print("# tileid ra dec", file=fx)
            for i in range(ntiles):
                print(f"{tileid[i]} {ra[i]:.6f} {dec[i]:.6f}", file=fx)
        files.append(output + '.dat')
        if verbose:
            print(f"Wrote {output}.dat")

    return dict(ntiles=ntiles, tilearea=actual_tilearea, output=output,
                tileid=tileid, ra=ra, dec=dec, files=files)


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate a high-density all-sky grid using a Fibonacci Sphere Lattice."
    )
    parser.add_argument("-o", "--output", type=str, help="output basename")
    parser.add_argument("-a", "--tilearea", type=float, required=True,
                        help="Desired average area coverage per point in square degrees (e.g., 0.25)")
    parser.add_argument("--no-ascii", action="store_true", help="Skip writing the .dat file")
    return parser.parse_args(argv)


def main(argv=None):
    opts = parse_arguments(argv)
    try:
        fibonacci_tiling(opts.tilearea, output=opts.output,
                         write_ascii=not opts.no_ascii)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
