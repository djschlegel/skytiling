#!/usr/bin/env python3
"""
Build the initial all-sky Fibonacci tiling for a given camera.

Given the camera astrometric file (one image HDU per CCD, each with a WCS
describing where that CCD lands on the sky for a reference pointing) and the
desired mean number of passes over the sky, this computes the on-sky area of
the camera footprint, sets

    tile_area = camera_area / npasses

and calls the Fibonacci lattice generator with that tile area.  The resulting
``allsky-fib-<ntiles>-<area>.fits`` is the input for ``optimize_tiling.py``.

Examples
--------
CFHT MegaCam (1.002 deg^2 footprint), four passes::

    make_initial_tiling.py -c data/cfht/cfht_40ccds.fits -n 4

which gives 0.2505 deg^2 per tile and 164,681 tiles.  The CFHT solution in
``data/cfht/`` was generated with ``--tilearea 0.25`` (165,012 tiles), which
is ``-n 4.008``; use ``--tilearea`` to request an exact area instead.
"""
import argparse
import os
import sys

# Allow running this file directly from a source checkout
if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from skytiling.fibonacci_tile import fibonacci_tiling, SPHERE_AREA_DEG2
from skytiling.optimize_tiling import CameraGeometry


def camera_area(camera_fits, verbose=False):
    """On-sky area of the camera footprint [deg^2] (sum over CCDs)."""
    cam = CameraGeometry(camera_fits, verbose=verbose)
    return cam.total_area_deg2, cam


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the initial Fibonacci tiling for a camera and a "
                    "mean number of passes.")
    parser.add_argument("-c", "--camera", required=True,
                        help="Camera astrometric FITS file (one image HDU per CCD with WCS)")
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("-n", "--npasses", type=float,
                   help="Mean number of passes over the sky (tile area = camera area / npasses)")
    g.add_argument("-a", "--tilearea", type=float,
                   help="Request an exact tile area [deg^2] instead of deriving it from npasses")
    parser.add_argument("-o", "--output",
                        help="Output basename (default allsky-fib-<ntiles>-<area>)")
    parser.add_argument("--no-ascii", action="store_true",
                        help="Only write the FITS file, not the .dat file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print the per-CCD footprint analysis")
    return parser.parse_args(argv)


def main(argv=None):
    opts = parse_arguments(argv)

    area, cam = camera_area(opts.camera, verbose=opts.verbose)
    print(f"Camera: {opts.camera}")
    print(f"  {len(cam.ccd_areas_deg2)} CCDs, footprint area = {area:.4f} deg^2, "
          f"max radius = {cam.max_radius * 180 / 3.141592653589793:.4f} deg")

    if opts.npasses is not None:
        if opts.npasses <= 0:
            sys.exit("ERROR: --npasses must be positive")
        tilearea = area / opts.npasses
        print(f"  Requested mean passes = {opts.npasses:g}  ->  tile area = {tilearea:.5f} deg^2")
    else:
        tilearea = opts.tilearea
        print(f"  Requested tile area = {tilearea:.5f} deg^2  ->  mean passes = {area / tilearea:.4f}")

    result = fibonacci_tiling(tilearea, output=opts.output, write_ascii=not opts.no_ascii)

    print(f"  {result['ntiles']} tiles, actual tile area = {result['tilearea']:.5f} deg^2, "
          f"actual mean passes = {area / result['tilearea']:.4f}")
    print(f"Next step: optimize_tiling.py -c {opts.camera} -t {result['output']}.fits ...")


if __name__ == "__main__":
    main()
