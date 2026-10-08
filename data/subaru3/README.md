# Subaru HSC + MBQ1, measured 8.18' cross model (rebuilt flats)

Third throughput model for HSC with the MBQ1 quadrant filter.  `data/subaru/`
was built from dome flats whose `NO_DATA` mask took the filter-holder cross to
be 10.4' wide and had no centre CCDs; `data/subaru2/` kept those flats and
*modelled* the band out to an 8.18' cross as a vignetting ramp.  This model
replaces the modelling with a **measurement**: Hironao Miyatake rebuilt the
flats with an 8.18' cross mask (half-width 4.09') and with the centre CCDs
included, and regenerated the throughput file with the same
`flats_to_throughput` defaults as before (0.1' cells, ±50', normalised at
r < 35', illumination trend fitted at 12' < r < 38').  The file is his
`hsc_mbq1_throughput_8p18.fits` from `hsc-niji-survey/wide-observing-strategy`
(`wide_dithering/`), 2026-10-08, unchanged; it has the same HDUs and header
keywords as the earlier files, so `skytiling` reads it without changes.

## Files

| File | Description |
|---|---|
| `hsc_mbq1_throughput.fits` | Camera + throughput model (`DETECTORS`, `FP2FA`, `DETFLATS` (94 rows), `GOOD_<f>`, `THRU_<f>`, `RADIAL`). |
| `hsc_mbq1_throughput_<f>.png` | One page per sub-filter from `plot_throughput --compare ../subaru2/hsc_mbq1_throughput.fits`: the filter's quadrant at full resolution, the throughput profiles toward the cross (solid; `subaru2` dotted), the cell-by-cell difference from `subaru2` (green outline = newly illuminated cells), and the radial profiles. |
| `anneal/` | Tilings annealed against this model (see `anneal/README.md`). |

The pointing list is still `../subaru/hsc_niji_wide_pointings.csv`.

## What changed relative to `subaru2`

* **Centre CCDs are back.** 413: `1_08`, 439: `1_16`, 465: `0_08`, 490: `0_16`
  — the CCD of each quadrant that touches the field centre — are now
  illuminated inward to r ≈ 5.9' (before: nothing inside r < 10').  They add
  0.011-0.012 deg² per sub-filter.  Their flats contain a light leak (a bright
  arc around the cross corner at r = 6-12', up to 3.4x the surroundings, in
  413/439/465; 490 is clean) that Hironao removed by fitting a smooth quadratic
  in (|x|, |y|) per CCD and amplifier at r > 13.5' and extrapolating inward, so
  **the throughput of the three centre CCDs at r < 13.5' is modelled, not
  measured** (extrapolation bias < 1%, scatter 1-4% where it could be tested).
  Inside r < 12' the radial illumination trend is likewise an extrapolation of
  the 12-38' fit, so the throughput of 0.78-0.84 at r ≈ 6.5' is a mix of the
  real cross penumbra and the assumed trend.  The leak itself was expected to
  be fixed at the telescope from 2026-10-10.
* **The cross band is measured.** The nearest illuminated cell is now at
  |x|, |y| = 4.05' (`subaru2`: 4.15' by construction).  The measured shadow is
  a soft penumbra rather than a step: in a raw MBQ1 flat the transmission
  across the arms is 0.98 / 0.95 / 0.90 / 0.75 / 0.50 at full widths of
  10.7 / 10.2 / 9.0 / 7.2 / 4.9' (vertical arm) and 16.3 / 12.3 / 9.7 / 6.3 /
  3.5' (horizontal arm) — equivalent opaque widths 5.1' and 4.3'.  So the
  10.4' mask of `subaru` was where the dimming is only 2-5%, and 8.18' is
  where it is 15-20%.  Compared with the `subaru2` ramp (0.8 at the opaque
  edge), the measured band is a little darker along the horizontal arm at
  large |x| and a little brighter near the centre (see the difference panels).
* `1_34` and `1_47` are still missing from the 413 set; dead amplifiers, the
  field edge and the CCD gaps are as before.

| | 413 | 439 | 465 | 490 |
|---|---|---|---|---|
| illuminated / throughput-weighted area, `subaru` | 0.2703 / 0.2583 deg² | 0.3020 / 0.2967 | 0.3005 / 0.2936 | 0.2903 / 0.2835 |
| `subaru2` | 0.2901 / 0.2748 | 0.3219 / 0.3141 | 0.3204 / 0.3111 | 0.3102 / 0.3013 |
| `subaru3` | 0.3012 / 0.2851 | 0.3328 / 0.3256 | 0.3302 / 0.3205 | 0.3188 / 0.3104 |

(illuminated = cells with `GOOD` > 0.5; the earlier READMEs counted `GOOD` > 0.99, hence slightly different numbers there.)

## Effect on the existing tilings

Evaluated on 8M independent randoms, footprint interior; per sub-filter:
mean, rms (rms/mean), fraction < 1.5, fraction < 0.5 exposures.  No re-annealing.

| | 413 | 439 | 465 | 490 |
|---|---|---|---|---|
| Atsushi phase 1 | 3.12, 1.12 (0.36), 6.9%, 0.45% | 3.54, 1.04 (0.29), 2.6%, 0.16% | 3.51, 1.06 (0.30), 2.8%, 0.36% | 3.39, 1.00 (0.30), 2.9%, 0.26% |
| `subaru2` adopted tiling (`../subaru2/anneal/nersc_atsushi2`) | 2.84, 0.91 (0.32), 6.5%, 0.49% | 3.24, 0.93 (0.29), 2.9%, 0.15% | 3.19, 0.89 (0.28), 2.9%, 0.14% | 3.09, 0.90 (0.29), 3.7%, 0.20% |

The centre CCDs help Atsushi's regular pattern a great deal — his four dithers
are sized to fill the cross, and the new area sits exactly there — so the
two tilings are now much closer than under the earlier models (rms/mean
0.36 / 0.29 / 0.30 / 0.30 vs 0.32 / 0.29 / 0.28 / 0.29; the annealed tiling
still has lower absolute rms in every sub-filter, but its mean inside the
interior is 8% lower because its tiles were pulled outward to cover the
footprint edges).  The tiling re-annealed against this model (`anneal/nersc_atsushi3/`,
`nersc/anneal_hsc_resume.sh`, 200 iterations) gains a further 0.7%: interior rms/mean
0.32 / 0.28 / 0.28 / 0.29; see `docs/hsc_mbq1_anneal_summary_subaru3.md`.

## Rectangular footprint: the DESI Run 2 Niji stripes

The DESI Run 2 imaging plan divides the survey into declination stripes 2.8
deg apart (`../desi_run2_stripes.csv`, its Table 2).  Four are to be imaged by
HSC-Niji + VST: NGC-5 and NGC-6 (adjacent, a 5.6-deg band at Dec -2.12 to
3.48) and SGC-1 and SGC-2 (adjacent, Dec -0.91 to 4.69).  `niji_stripes_2p8.csv`
gives their unique 2.8-deg bands with the tabulated RA bounds (1134 deg² in
all) — the region where the medium-band coverage should be uniform, with the
hard Dec edges abutting the stripes imaged by other telescopes.  (The plan's
tabulated Dec bounds are 3.5 deg tall, i.e. 0.35 deg of padding on each side;
pass those instead if the padding is to be covered at full depth too.)

`optimize_tiling_maps --rects niji_stripes_2p8.csv` uses this union of
rectangles as the footprint (randoms uniform inside it; "interior" statistics
exclude points within `--interior-margin` 0.75 deg of its boundary).
`--init-lattice 1.299,1.125` seeds pointing centres on rows of constant Dec
with the HSC-SSP Wide spacing, fitted to each band (5 rows per 5.6-deg band,
1.12 deg apart, with the number of columns rounded so the lattice is centred
in each row's RA extent), and `--init-offsets` adds Atsushi's four dithers:
1.45 deg² per centre as in his layout (784 centres / 3136 tiles on the DESI bounds,
637 / 2548 on the trimmed stripes).  Alternatively `--init-centers
../subaru/hsc_niji_wide_pointings.csv --margin 0` starts from the HSC-SSP Wide
pointings themselves, keeping the 640 whose centres fall inside the trimmed stripes
(five rows per band in each cap, e.g. Dec -1.3, -0.18, 0.95, 2.08, 3.2 in the NGC): the
two starts have the same density and, with the four dithers, the same coverage
statistics to within 1-2% (lattice: 6.1% of the coverage outside the footprint, inside
rms/mean 0.39 / 0.34 / 0.35 / 0.34; SSP pointings: 6.7%, 0.40 / 0.34 / 0.35 / 0.35 — the
SSP rows sit 0.28 deg below the top of the NGC band and 0.82 deg above its bottom, the
lattice rows 0.56 deg from both edges).  The output tables carry a `pid` column (index
of the undithered centre; with `--init-centers`, the row index in that CSV) and, when
the CSV has names, a `pname` column (e.g. `W_010001`), so the four dithers of a
pointing can be scheduled together and traced to the SSP pointing.

To stop the anneal pushing coverage outside the footprint, `--objective cv
--cv-power p` minimises Σ_f rms_f / mean_f^p instead of Σ_f N Var(c_f).
Because the total coverage of a fixed set of exposures is conserved, the mean
inside the footprint rises only when coverage is pulled back in from outside,
so the cv objective charges for spill; p sets how strongly.  On a 15 x 5.6 deg
test rectangle (232 tiles, 25 local iterations from the lattice start):

| objective | mean inside (all randoms) | rms/mean, all randoms | interior mean | interior rms/mean | interior < 1.5 exp. |
|---|---|---|---|---|---|
| start (lattice + dithers) | 2.96 / 3.37 / 3.32 / 3.21 | 0.41 / 0.36 / 0.37 / 0.36 | 3.25 / 3.73 / 3.67 / 3.55 | 0.35 / 0.29 / 0.30 / 0.29 | 5.9 / 2.4 / 2.8 / 2.9% |
| var (p = 0) | 2.80 / 3.19 / 3.14 / 3.04 | 0.33 / 0.30 / 0.30 / 0.31 | 2.90 / 3.32 / 3.26 / 3.17 | 0.32 / 0.28 / 0.28 / 0.29 | 5.9 / 2.6 / 2.8 / 3.5% |
| cv, p = 1 | 2.92 / 3.32 / 3.28 / 3.17 | 0.33 / 0.30 / 0.30 / 0.31 | 3.07 / 3.51 / 3.46 / 3.35 | 0.31 / 0.28 / 0.28 / 0.29 | 4.7 / 2.1 / 2.0 / 2.7% |
| cv, p = 2 | 3.02 / 3.44 / 3.38 / 3.28 | 0.34 / 0.31 / 0.31 / 0.32 | 3.21 / 3.68 / 3.62 / 3.51 | 0.31 / 0.28 / 0.28 / 0.29 | 3.5 / 1.6 / 1.7 / 2.2% |

The variance objective lost 5% of the mean to spill; p = 1 reaches the same
uniformity with the spill mostly held, and p = 2 pulls the tiles in further
(interior depth +10% over var, fewer holes) at a small cost in the global
rms/mean, which now includes a steeper edge roll-off.  `--keep-inside d`
additionally rejects moves that take a tile centre more than d deg beyond the
footprint (not needed in this test: no centre got that far).
`nersc/anneal_hsc_stripes.sh [p]` runs the full four-stripe problem on the trimmed footprint.

### How the stripes compare with the HSC-Niji wide footprint

`niji_stripes_vs_hsc.png` overlays the stripes on the hex cells of the 802
HSC-Niji wide pointings (HSC-SSP Wide DR4 full-colour pointings, Atsushi's
footprint).  Of the stripes' 1135 deg², 905 lie inside those cells and 230 outside;
266 deg² of the HSC footprint (its top row in the NGC at Dec 3.5-5.1, and the SGC
area at Dec > 4.7 and the block at RA 30-39, Dec -2 to -6) is not in any stripe.
Per stripe (deg², and the RA / Dec extent of the HSC cells within the stripe's Dec band):

| stripe | area | inside HSC | outside HSC | HSC RA extent | HSC Dec extent |
|---|---|---|---|---|---|
| NGC-5 | 345 | 273 | 72 | 128.0-226.0 | 0.72-3.47 |
| NGC-6 | 280 | 250 | 31 | 128.2-226.0 | -2.03-0.67 |
| SGC-1 | 248 | 191 | 57 | 330.2-39.7 | 1.92-4.67 |
| SGC-2 | 261 | 192 | 70 | 330.2-39.7 | -0.88-1.87 |

The excess is almost all in RA: NGC-5 runs 24 deg past the HSC data at the high-RA end,
NGC-6 2 deg, and the SGC stripes extend 16 deg below and 5-8 deg above it.  In Dec the
NGC band fits the HSC rows well (its bottom edge at -2.12 runs along the zig-zag bottom of
the lowest row of cells, so a 0.4-deg strip there is only half covered); the SGC band
sits inside the HSC Dec range.  **`niji_stripes_2p8_hsc.csv` — the adopted footprint —** is the same four stripes with
the RA bounds trimmed to the HSC extent (NGC 128.0-226.0, SGC 330.2-39.7): 937 deg², of
which 905 inside HSC and 32 outside (the NGC-6 bottom strip); 637 lattice centres, 2548
tiles with Atsushi's four dithers.  `niji_stripes_hsc_vs_hsc.png` is the same overlay for the
trimmed stripes.

| stripe (trimmed) | RA | Dec (2.8-deg band) | area deg² |
|---|---|---|---|
| NGC-5 | 128.0-226.0 | 0.68-3.48 | 274.2 |
| NGC-6 | 128.0-226.0 | -2.12-0.68 | 274.4 |
| SGC-1 | 330.2-39.7 | 1.89-4.69 | 194.3 |
| SGC-2 | 330.2-39.7 | -0.91-1.89 | 194.6 |
