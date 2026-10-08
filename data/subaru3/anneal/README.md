# Annealed HSC MBQ1 tiling for the measured 8.18' cross model

Output of `optimize_tiling_maps` against `../hsc_mbq1_throughput.fits` (see
`../README.md`).  File format as in `../../subaru/anneal/README.md`: one row
per pointing centre, observed at rotation offsets `rot`, `rot+90`, `rot+180`,
`rot+270`.

| directory | contents |
|---|---|
| `nersc_atsushi3/` | Perlmutter run of `nersc/anneal_hsc_resume.sh` (200 iterations) starting from the `data/subaru2` solution (`../../subaru2/anneal/nersc_atsushi2/tiles_0200.fits`, itself annealed from the `data/subaru` solution, itself from Atsushi's phase-1 layout; 3208 pointings): `tiles_0000.fits` (start), `tiles_0200.fits` (final, full precision), `convergence.txt`, and the adopted `hsc_mbq1_tiling_subaru3_final.{fits,csv}` (`tileid, ra, dec, rot`; RA/Dec to 1e-3 deg, ROT to 0.1 deg). *To be filled in when the run completes.* |

Run: 4M randoms (`hsc_fib/fib_randoms_4000000.fits`, seed 1), `--max-drift 0.75
--moves-per-tile 4 --delta 0.05 --delta-rot 5 --shrink 0.99 --seed 3`, 64 workers;
the drift cap is still measured from Atsushi's original positions (`ra0, dec0`).
