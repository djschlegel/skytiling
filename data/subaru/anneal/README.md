# Annealed HSC MBQ1 tilings of the HSC-Niji wide footprint

Outputs of `optimize_tiling_maps` (see `../README.md` and
`../../../docs/hsc_mbq1_anneal_summary.md`).  Every tiles file has one row per
pointing centre with columns `tileid, ra, dec, rot, ra0, dec0`: the pointing is
observed at instrument rotation offsets `rot`, `rot+90`, `rot+180`, `rot+270`
(degrees; `ROTS` header keyword), and `ra0, dec0` is the starting position the
drift was limited from.  Header keywords `MEAN_<f>` / `RMS_<f>` give the
coverage statistics over all randoms at that iteration.

| directory | contents |
|---|---|
| `nersc_atsushi/` | Perlmutter run starting from Atsushi's phase-1 dither pattern (4 dithers on each of the 802 wide pointings = 3208 tiles): `tiles_0000.fits` (the start, rot = 0) and `tiles_0300.fits` (**final**), plus `convergence.txt` (per-iteration mean and rms per filter). |
| `nersc_fibonacci/` | Perlmutter run starting from a Fibonacci lattice of 3207 tiles inside the footprint: `tiles_0000.fits` (start) and `tiles_0276.fits` (**final**), `convergence.txt`. |
| `local_40iter/` | Earlier 40-iteration runs on a laptop (2M randoms): from Atsushi's layout (`atsushi4_*`), from a Fibonacci start without (`fibonacci_*`) and with (`fibonacci_rot_*`) rotation annealing; logs and coverage plots included. |

Both NERSC runs used 4M randoms (`fib_randoms_4000000.fits`, 64 MB, not in
git; regenerate with `--num-randoms 4000000 --seed 1`), `--max-drift 0.75
--moves-per-tile 4 --delta 0.1 --delta-rot 10 --shrink 0.99`, 64 workers.
