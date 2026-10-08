# Annealed tiling for HSC + MBQ1 on the HSC-Niji wide footprint

*D. Schlegel, 2026-10-08. Code: `skytiling` (`optimize_tiling_maps`), data in `data/subaru/`.*

**Setup.** The camera is described by per-sub-filter maps of relative throughput on the focal plane,
built from Hironao's MBQ1 dome flats (gain-corrected, sky-area-per-pixel and dome-illumination
gradient divided out, ≈1 over the bulk of the field, 0 where the `NO_DATA` mask flags the
~10.4′-wide filter-holder cross, the vignetted edge beyond 46.5′, dead amplifiers and the two CCDs
absent from the 413 flats). Every pointing is observed at four instrument rotations 90° apart, and
the coverage of a sky point in sub-filter *f* is the sum of throughput over all exposures. Starting
from a Fibonacci lattice of **3207 pointings** inside the 1194 deg² footprint (0.75° discs around the
802 HSC-SSP Wide centres at −6° < Dec < 6.4°) — the same number of exposures as Hironao's phase-1
pattern (4 dithers × 802 pointings × 4 PAs) — simulated annealing moves each pointing centre (≤ 0.75°)
and the rotation offset of its 4-PA set to minimise the summed per-filter coverage variance over 4M
random points. 277 iterations on one Perlmutter node (~6 s per iteration with 64 workers; step size
0.1° → 0.03°, 10° → 3° at iteration 77).

**Result** (footprint interior, 68% of the area, away from the stripe edges; coverage = effective
number of full-throughput exposures):

| | 413 | 439 | 465 | 490 |
|---|---|---|---|---|
| Hironao phase 1: mean ± rms (rms/mean) | 2.85 ± 1.19 (0.42) | 3.24 ± 1.12 (0.34) | 3.24 ± 1.16 (0.36) | 3.11 ± 1.10 (0.35) |
| &nbsp;&nbsp;&nbsp; fraction < 1.5 / < 0.5 | 13.0% / 2.1% | 6.5% / 0.8% | 7.4% / 1.2% | 7.7% / 0.9% |
| Fibonacci lattice, un-annealed | 2.79 ± 1.08 (0.39) | 3.19 ± 1.08 (0.34) | 3.17 ± 1.15 (0.36) | 3.06 ± 1.04 (0.34) |
| **Annealed (positions + rotations)** | **2.56 ± 0.86 (0.34)** | **2.94 ± 0.88 (0.30)** | **2.91 ± 0.85 (0.29)** | **2.81 ± 0.84 (0.30)** |
| &nbsp;&nbsp;&nbsp; fraction < 1.5 / < 0.5 | 10.2% / 0.8% | 4.9% / 0.3% | 4.7% / 0.2% | 5.8% / 0.3% |

For the same number of exposures the annealed tiling lowers the rms of the coverage by 24–27% in
every filter (rms/mean by 13–19%), cuts the deep holes (coverage < 0.5) by a factor 2.8–5, and
removes the periodic pattern of Hironao's layout (Figure 1). The mean is ~10% lower because pointings
near the stripe boundary move outward to flatten the edge and the coverage they put outside the
footprint is not counted; in the interior the annealed solution is more uniform at every quantile.
413 remains the worst filter in both cases because of its two missing CCDs. The 0/1 (`GOOD`) and
throughput-weighted statistics differ by only 2–4%.

**Rotations.** Allowing the 4-PA set of each pointing to turn is worth a few per cent in rms on top
of the position moves. The annealed offsets (Figure 2, lower left) have an rms of 17°, with 57% of
pointings turned by more than 10° and 9% by more than 30°; the distribution is broad and single-peaked
with a small positive mean (+3.5°), i.e. the solution does not favour a particular grid of angles. The
pointings themselves moved by a median 0.19° from the lattice (max 0.75° = the cap, reached by only 3
pointings), so the Fibonacci lattice is a good starting point and the result is set by the camera
masks and the tile density rather than by the start.

**Convergence.** The objective fell from 1.58 to 1.24 (sum over filters of rms/mean); 85% of the gain
came in the first 40 iterations, and iterations 207–276 gained only 0.005 (Figure 2, upper left). At the
final slope of 6e-5 per iteration, another 300 iterations would buy well under 1%, so the run is
converged for practical purposes at this tile density. Larger gains would come from more tiles (lower rms/mean at
higher mean coverage) or from a different objective (e.g. penalising holes rather than variance).

![](figures/hsc_fib_coverage_compare.png)

*Figure 1. Throughput-weighted coverage per sub-filter in a 4° × 4° field around (RA, Dec) = (180°, 0°)
for Hironao's phase-1 pattern (top) and the annealed tiling (bottom, iteration 276), same colour scale.*

<img src="figures/hsc_fib_anneal_stats.png" style="width:84%">

*Figure 2. Convergence of the objective; distribution of MBQ1-413 coverage over the footprint interior;
annealed rotation offsets; and displacement of each pointing from its initial lattice position.*

Files: `hsc_fib/fib_tiles_0276.fits` (columns `ra, dec, rot, ra0, dec0`; each row is one pointing to be
observed at `rot`, `rot+90`, `rot+180`, `rot+270`), `hsc_fib/fib_randoms_4000000.fits`. The mapping
of `rot` onto `INSROT_PA` (sign and zero point) still has to be calibrated against the instrument.
