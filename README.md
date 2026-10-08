# skytiling

Uniform all-sky survey tilings for arbitrary camera footprints.

Given a camera's astrometric layout (one WCS per CCD), `skytiling` produces a
set of pointing centres that covers the whole sky with a chosen mean number of
passes and as uniform a coverage as possible.  It was developed for the CFHT
MegaCam 40-CCD mosaic, and is written so the same code runs unchanged for other
cameras (Subaru HSC, Rubin LSSTCam, VST OmegaCAM, ...) given their camera file.

## How it works

1. **Camera footprint** — `cfht_40ccds.fits`-style file: one image HDU per
   CCD whose header carries a WCS (TAN-SIP or similar) for a reference
   pointing.  Only the headers and pixel dimensions are used; the pixel data
   can be anything.  The on-sky area of each CCD is computed from the WCS
   corners (spherical quads).
2. **Initial tiling** — `make_initial_tiling` divides the footprint area by the
   requested mean number of passes to get a tile area, then lays down a
   Fibonacci sphere lattice with that density (`fibonacci_tile`).  The
   Fibonacci lattice is equal-area and very uniform, so it is already a good
   starting point.
3. **Optimisation** — `optimize_tiling` drops a catalogue of random points on
   the sky, counts how many times each is covered by a CCD, and runs
   simulated annealing on the tile centres to minimise the RMS of the coverage
   counts (i.e. make the number of passes as uniform as possible).  It
   checkpoints every iteration (`output_tiles_NNNN.fits`) and can resume from
   any checkpoint.

## Installation

```
git clone https://github.com/<you>/skytiling
cd skytiling
pip install -e .          # installs the three console scripts
```

or run the scripts straight from the checkout with `bin/<script>` (no
install needed; requires numpy, scipy, astropy, matplotlib).

## Usage

### 1. Initial tiling from a camera file and a mean number of passes

```
make_initial_tiling -c data/cfht/cfht_40ccds.fits -n 4
```

prints the footprint area and writes `allsky-fib-<ntiles>-<area>.fits`
(and `.dat`).  Use `-a 0.25` instead of `-n` to request an exact tile area,
and `-v` to see the per-CCD areas.  The lower-level generator is also
available directly as `fibonacci_tile -a <deg2>`.

### 2. Optimise the tiling

```
optimize_tiling \
    -c data/cfht/cfht_40ccds.fits \
    -t allsky-fib-165012-0.250.fits \
    --num-randoms 10000000 --max-drift 1.5 --delta 0.1 --temp 0.0 \
    --iters 100 --moves-per-tile 4 \
    --plot --ra-center 10 --dec-center 0 --diameter 5 \
    --write-coverage
```

Each iteration writes `output_tiles_NNNN.fits` (tile table, plus an optional
coverage HDU with `--write-coverage`) and, with `--plot`,
`output_coverage_NNNN.png` / `output_layers_NNNN.png`.  Pass a checkpoint as
`-t output_tiles_0058.fits` to resume; the iteration counter continues from
the file name.  `optimize_tiling --help` lists all options, including
`--fixed-{ra,dec}-{min,max}` to freeze tiles outside a region and
`--workers N` for multiprocessing (on macOS use `--force-fork`).

## Layout

```
skytiling/            Python package
  fibonacci_tile.py     Fibonacci-lattice all-sky point generator
  make_initial_tiling.py  camera file + mean passes -> initial tiling
  optimize_tiling.py    simulated-annealing coverage optimiser
  lsst_camera.py        read the camera geometry embedded in LSST-pipeline (HSC, LSSTCam) files
  flats_to_throughput.py  LSST-pipeline flats -> field-angle throughput/footprint maps
  eval_dither_tiling.py   coverage statistics of a hex-grid + dither + rotation tiling with those maps
  optimize_tiling_maps.py annealing optimiser for throughput-map cameras (multi-filter, rotations, footprint)
bin/                  wrappers to run the above from a source checkout
data/cfht/            CFHT MegaCam camera file, initial tiling, final solution
data/subaru/          Subaru HSC + MBQ1 quadrant filter: camera and throughput model
docs/                 run notes
nersc/                Perlmutter batch scripts and notes for the HSC runs
```

Each `data/<camera>/` directory has its own README describing the files.
Large or reproducible files (full-image camera files, random catalogues,
intermediate iterations) live in `data/<camera>/scratch/`, which is
git-ignored.

## Camera file format

A multi-extension FITS file with one image HDU per CCD.  Each HDU needs
`NAXIS1`/`NAXIS2` (the CCD size in pixels) and a WCS that
`astropy.wcs.WCS` can read; `CRVAL1`/`CRVAL2` of the first CCD define the
reference boresight.  Pixel values are ignored, so a compact file with
tile-compressed constant images (as in `data/cfht/cfht_40ccds.fits`) is
equivalent to a full image mosaic.

## License

BSD 3-Clause; see `LICENSE`.
