# Subaru HSC tiling

Inputs in preparation.  Expected contents, following `data/cfht/`:

- `subaru_<n>ccds.fits` — camera astrometric file, one image HDU per CCD with a WCS
  for a reference pointing (`CRVAL1/2` of the first HDU = boresight).
- `allsky-fib-<ntiles>-<area>.fits` — initial tiling from
  `make_initial_tiling -c subaru_<n>ccds.fits -n <passes>`.
- `output_tiles_NNNN.fits` — final optimized solution.
- `scratch/` (git-ignored) — randoms, intermediate iterations, full-image files.
