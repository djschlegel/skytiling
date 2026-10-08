# Annealed HSC MBQ1 tiling for the 8.18' cross model

Output of `optimize_tiling_maps` against `../hsc_mbq1_throughput.fits` (see
`../README.md` and `../../../docs/hsc_mbq1_anneal_summary_subaru2.md`).
File format as in `../../subaru/anneal/README.md`: one row per pointing
centre, observed at rotation offsets `rot`, `rot+90`, `rot+180`, `rot+270`.

| directory | contents |
|---|---|
| `nersc_atsushi2/` | Perlmutter run of 200 iterations starting from the `data/subaru` solution (`../../subaru/anneal/nersc_atsushi/tiles_0300.fits`, itself annealed from Atsushi's phase-1 layout; 3208 pointings): `tiles_0000.fits` (start), `tiles_0200.fits` (final, full precision), `convergence.txt`. **`hsc_mbq1_tiling_subaru2_final.{fits,csv}` is the adopted solution** (`tileid, ra, dec, rot`; RA/Dec to 1e-3 deg, ROT to 0.1 deg). |

Run: 4M randoms (`hsc_fib/fib_randoms_4000000.fits`, seed 1), `--max-drift 0.75
--moves-per-tile 4 --delta 0.05 --delta-rot 5 --shrink 0.99 --seed 2`, 64 workers.
