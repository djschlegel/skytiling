# Subaru HSC + MBQ1, narrower cross (8.18')

Same camera and throughput model as `data/subaru/`, with one change: the
opaque part of the filter-holder cross is taken to be **8.18' wide**
(|x| or |y| < 4.09') instead of the 10.4' that the `NO_DATA` mask of the dome
flats implies.  The band between 4.09' and 5.2' from each axis — real CCD
area that the flats mask — is treated as illuminated but vignetted: the
throughput ramps linearly, along each row/column, from the flat's own value
just outside the old mask edge down to **0.8 at the opaque edge**.  The flats
already show this ramp beginning further out (in 413 the throughput falls from
1.0 at ~11' from an axis to 0.86 at 5.2', 0.82 next to the x axis; in 439/465
from 1.0 at ~8.5' to ~0.96; 490 shows almost no ramp), so the band continues
each quadrant's measured trend.  Everything else (field edge, dead
amplifiers, missing 413 CCDs, CCD gaps, radial throughput) is unchanged.

```
adjust_cross -i ../subaru/hsc_mbq1_throughput.fits -o hsc_mbq1_throughput.fits \
    --old-half-width 5.2 --new-half-width 4.09 --band-throughput 0.8 --band-profile ramp --plot
```

`hsc_mbq1_throughput.png` shows the central 32' of each sub-filter at full
map resolution (0.1') with the band outlined in red, and the throughput
profiles toward the cross along x and y.  The quadrant-to-quadrant difference
in how early the vignetting starts (strong in 413, absent in 490) comes
straight from the flats and may be worth checking against the filter-holder
geometry.

| | 413 | 439 | 465 | 490 |
|---|---|---|---|---|
| illuminated area per exposure, `subaru` | 0.2705 deg² | 0.3020 | 0.3008 | 0.2906 |
| illuminated area, `subaru2` | 0.2900 | 0.3215 | 0.3203 | 0.3101 |
| throughput-weighted area, `subaru2` | 0.2770 | 0.3165 | 0.3134 | 0.3036 |

The band adds 0.020 deg² (7%) of illuminated area per sub-filter.  `hsc_mbq1_throughput.png`
shows the maps with the band outlined in red.  The footprint pointing list is the same
`../subaru/hsc_niji_wide_pointings.csv`.

## Effect on the existing tilings

Evaluated on 4M independent randoms (footprint interior; mean, rms (rms/mean), fraction
< 1.5, fraction < 0.5), without re-annealing:

| | model | 413 | 439 | 465 | 490 |
|---|---|---|---|---|---|
| Atsushi phase 1 | subaru | 2.82, 1.19 (0.42), 13.0%, 2.1% | 3.22, 1.12 (0.35), 6.5%, 0.8% | 3.21, 1.17 (0.36), 7.5%, 1.2% | 3.08, 1.10 (0.36), 7.7%, 1.0% |
| Atsushi phase 1 | subaru2 | 3.00, 1.15 (0.38), 9.1%, 0.8% | 3.40, 1.09 (0.32), 4.0%, 0.3% | 3.40, 1.11 (0.33), 4.2%, 0.5% | 3.27, 1.04 (0.32), 4.1%, 0.3% |
| annealed (`subaru/anneal/nersc_atsushi`) | subaru | 2.58, 0.90 (0.35), 10.5%, 0.9% | 2.96, 0.90 (0.31), 5.0%, 0.3% | 2.93, 0.87 (0.30), 4.9%, 0.3% | 2.83, 0.87 (0.31), 6.0%, 0.4% |
| annealed (`subaru/anneal/nersc_atsushi`) | subaru2 | 2.74, 0.91 (0.33), 7.9%, 0.6% | 3.13, 0.91 (0.29), 3.6%, 0.2% | 3.10, 0.89 (0.29), 3.5%, 0.2% | 3.00, 0.88 (0.30), 4.3%, 0.2% |

The narrower cross helps every tiling (the mean rises 5-6% and the holes
shrink), and it helps Atsushi's regular pattern more than the annealed one —
the cross band is exactly where his dithers leave gaps — but the annealed
tiling remains ahead: rms/mean 0.33 / 0.29 / 0.29 / 0.29 vs 0.38 / 0.32 /
0.33 / 0.32, with ~15-25% fewer holes.  A tiling annealed against the
`subaru2` maps themselves (`optimize_tiling_maps -t data/subaru2/hsc_mbq1_throughput.fits ...`)
should do a little better still.
