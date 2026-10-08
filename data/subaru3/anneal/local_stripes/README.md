# Stripe-footprint anneals (local runs, 100 iterations)

Three anneals of the HSC MBQ1 tiling over the adopted stripe footprint
(`../../niji_stripes_2p8_hsc.csv`, 999 deg², 6 rectangles) against the `data/subaru3`
model, run on a laptop with 2M randoms, 100 iterations, `--delta 0.1 --shrink 0.98
--moves-per-tile 2 --delta-rot 10 --max-drift 0.75` (the commands are in
`nersc/anneal_hsc_stripes.sh` apart from the iteration count and random number).
A 151-iteration run of the `ssp_p2` case on 1M randoms gave the same statistics to
within 0.005 in rms/mean, so 100 iterations are converged for this purpose.

| run | start | objective | tiles | files |
|---|---|---|---|---|
| `ssp_p1` | 680 HSC-SSP Wide pointings inside the stripes x Atsushi's 4 dithers | rms/mean (p = 1) | 2720 | `ssp_p1_tiles_{0000,0100}.fits`, `ssp_p1.out` |
| `ssp_p2` | same | rms/mean² (p = 2) | 2720 | `ssp_p2_tiles_0100.fits`, `ssp_p2.out` |
| `lattice_p1` | 1.299 x 1.125 deg lattice fitted to the bands x 4 dithers | rms/mean (p = 1) | 2696 | `lattice_p1_tiles_{0000,0100}.fits`, `lattice_p1.out` |

Tables carry `pid` (and `pname` for the SSP start) identifying the undithered pointing.

## Evaluation on 6M independent randoms (footprint + 1.2 deg margin)

"spill" is the fraction of the total coverage that lands outside the stripes;
"inside" is all randoms in the stripes, "interior" those > 0.75 deg from the boundary.
Per sub-filter 413 / 439 / 465 / 490.

| run | spill [%] | mean inside | rms/mean inside | rms/mean interior | interior < 1.5 exposures [%] |
|---|---|---|---|---|---|
| start (SSP + dithers) | 6.2 / 6.4 / 6.3 / 6.4 | 2.98 / 3.38 / 3.33 / 3.23 | 0.40 / 0.34 / 0.35 / 0.34 | 0.35 / 0.29 / 0.30 / 0.29 | 6.5 / 2.5 / 2.8 / 2.7 |
| `ssp_p1` | 7.6 / 7.9 / 7.8 / 7.8 | 2.94 / 3.32 / 3.28 / 3.18 | 0.32 / 0.29 / 0.29 / 0.30 | 0.31 / 0.28 / 0.27 / 0.29 | 5.2 / 2.4 / 2.4 / 3.2 |
| `ssp_p2` | 5.3 / 5.5 / 5.4 / 5.4 | 3.02 / 3.41 / 3.37 / 3.27 | 0.33 / 0.30 / 0.29 / 0.30 | 0.31 / 0.28 / 0.27 / 0.29 | 4.5 / 2.0 / 2.0 / 2.7 |
| `lattice_p1` | 7.0 / 7.2 / 7.1 / 7.1 | 2.93 / 3.32 / 3.28 / 3.17 | 0.32 / 0.29 / 0.29 / 0.30 | 0.31 / 0.28 / 0.27 / 0.29 | 5.2 / 2.4 / 2.4 / 3.1 |

All three anneals reach the same uniformity (rms/mean 0.31-0.32 / 0.28 / 0.27 / 0.29
in the interior, from 0.35 / 0.29 / 0.30 / 0.29); the start barely matters.  p = 2
keeps 1-2% more coverage inside the stripes than p = 1 and has fewer deep holes
(4.5 vs 5.2% of the interior below 1.5 exposures in 413), at the cost of a slightly
steeper edge; it is the better choice if the stripes are the whole imaging area.

`hsc_stripes_edge_profiles.png`: coverage as a function of Dec across the NGC and SGC
bands.  With hard band edges the coverage at the boundary is 50-80% of the plateau
and reaches 90% of it only 0.4-0.8 deg inside, for every objective: a camera 1.5 deg
across cannot stop at a line, so full depth to the band edge requires imaging ~0.4 deg
beyond it (the padded 3.5-deg bounds of the DESI stripe table).  The objective only
decides how the roll-off is shared between the edge zone and spill.
`hsc_stripes_stats.png`: convergence and the 413 coverage distribution.
