# CFHT MegaCam tiling

Inputs and final solution for the CFHT MegaCam 40-CCD mosaic
(footprint 1.0020 deg², max radius 0.699°).

| File | Description |
|---|---|
| `cfht_40ccds.fits` | Camera astrometric file: 40 image HDUs with TAN-SIP WCS. Compact copy — headers identical to the original mosaic, pixel data replaced by tile-compressed constant images (optimize_tiling only uses headers and dimensions; verified to give identical geometry). |
| `allsky-fib-165012-0.250.fits` | Initial Fibonacci tiling: 165,012 tiles at 0.250 deg²/tile = 4.008 mean passes. Reproduce with `make_initial_tiling -c cfht_40ccds.fits -a 0.25`. |
| `output_tiles_0141.fits` | **Final solution** after 141 optimizer iterations: `tileid, ra, dec` for 165,012 tiles. The 300M-row coverage HDU has been stripped (see below). |
| `output_coverage_0141.png` | Coverage-count histogram / diagnostics at iteration 141. |
| `output_layers_0141.png` | Zoomed layer plot (RA 10°, Dec 0°, 5° diameter) at iteration 141. |
| `output_camera.png` | Camera footprint plot. |
| `run_cfht.sh` | The optimizer invocations used to produce the solution. |

## scratch/ (git-ignored)

Not in the repository; kept locally only:

- `cfht_40ccds_fullimage.fits` — the original 1.5 GB camera file with full float32 images.
- `output_tiles_0141_with_coverage.fits` — final solution including the coverage HDU
  (one int32 count per random, 30M randoms × 10 layers = 300M rows, 1.2 GB).
- `output_tiles_0000.fits`, `output_*_0000.png` — iteration-0 checkpoint and plots.
- `generated_randoms_10000000.fits`, `generated_randoms_30000000.fits` — random catalogues.

## Run history

Iterations 0–58 used 10M randoms with 4 moves per tile (≈13 min/iteration
on a Mac); iterations 59–141 continued from `output_tiles_0058.fits` with 30M
randoms (≈25 min/iteration).  See `run_cfht.sh` and `docs/cfht_run_notes.txt`.
