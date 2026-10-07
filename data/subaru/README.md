# Subaru HSC with the MBQ1 quadrant filter

Inputs for tiling the sky with Hyper Suprime-Cam and the **MBQ1** four-quadrant
narrow-band filter (413 / 439 / 465 / 490 nm, one per quadrant).  Each pointing
is to be observed at four instrument rotations (0, 90, 180, 270 deg) so every
sub-filter sweeps every quadrant; the filter holder's cross vignettes the field
between quadrants, and the field edge is vignetted, so coverage has to be
treated as a throughput between 0 and 1 rather than 0/1.

## Files

| File | Description |
|---|---|
| `hsc_mbq1_throughput.fits` | **Compact camera + throughput model** built by `flats_to_throughput` from Hironao Miyatake's dome flats (below). 3.8 MB. |
| `hsc_mbq1_throughput.png` | The maps in that file. |
| `flat_run2.tar.gz` | *(git-ignored, 3.9 GB)* Hironao's `flat_mbq1_{413,439,465,490}_run2_mask_nocenter` LSST-pipeline flats, 90 files. |
| `scratch/` | *(git-ignored)* the untarred flats and a clone of `hsc-niji-survey/wide-observing-strategy`. |

### `hsc_mbq1_throughput.fits`

Everything is in **field-angle coordinates**: the gnomonic (tangent-plane) offset
(x, y) from the boresight in the camera's own frame (LSST `FIELD_ANGLE`), in
arcmin for the maps and degrees for the detector corners.  The sky position
for a pointing at (RA, Dec, PA) is the TAN projection about the boresight of the
field angle rotated by PA.

| HDU | Contents |
|---|---|
| `DETECTORS` | 112 HSC detectors (104 science + guide/focus): name (`1_00` ...), pixel size `nx, ny`, the exact affine `focal plane [mm] = A[:, :2] @ (px, py) + A[:, 2]` (0-based pixels), and the field-angle corners of each CCD. |
| `FP2FA` | The optical distortion: `field_angle [rad] = sum coef * x_mm^px * y_mm^py`, 55 terms per axis (9th order), straight from the camera model embedded in the flats. |
| `DETFLATS` | Per (sub-filter, detector): fraction of illuminated pixels and median flat value. |
| `GOOD_<f>` | 1000x1000 map, 0.1'/cell, +-50': fraction of pixels illuminated by sub-filter `f` (not flagged `NO_DATA`). 0/1 away from edges; NaN where there is no detector (CCD gaps, outside the mosaic). |
| `THRU_<f>` | Same grid: relative throughput, ~1 over the bulk of the field, 0 where masked, NaN where no detector. Built as gain-corrected flat / (sky area per pixel from the distortion model) x `GOOD`, normalised per sub-filter at r < 35', divided by a linear-in-radius trend fitted at 12' < r < 38' (the dome-flat illumination gradient) and capped at 1.0. Departures from 1 are the edge vignetting (r > 43'), CCD-to-CCD QE differences of a few %, dead amplifiers and the masked cross/edges. |
| `RADIAL` | Azimuthal medians vs field radius per sub-filter (1' bins): the gain-corrected flat over the Jacobian (`median_raw`), the sky area per pixel relative to on-axis (`jacobian`), the fitted illumination trend (`model`) and the final throughput (`median_thru`). |

Read it with astropy; `skytiling.lsst_camera.LsstCamera.from_hdus(h['DETECTORS'], h['FP2FA'])`
rebuilds the geometry.

## What the flats are, and what was found in them

The tarball holds one LSST-pipeline (afw) `flat` per (sub-filter, CCD):
`flat_HSC_MBQ1_MBQ1_<ccd>_u_miyatake_flat_mbq1_<nm>_run2_mask_nocenter_...fits`,
each a MaskedImage (IMAGE / MASK / VARIANCE) of a 2048 x 4176 CCD, built from
`DATA-TYP = DOMEFLAT` exposures on 2026-06-17 (PROP-ID o26134).  Only the 21-23
CCDs under each quadrant are present (90 of the 104 science CCDs).  Points:

* **Camera geometry comes for free.** Every file embeds the full HSC `Camera`
  (per-CCD pixel->focal-plane affines and the focal-plane->field-angle polynomial)
  as AST text; `skytiling/lsst_camera.py` parses it with no LSST-stack
  dependency.  Plate scale 11.198"/mm = 0.1680"/pixel on axis; the CCD mosaic
  reaches r = 47'.  There is no sky WCS in the flats (dome flats), and none is needed.
* **The quadrant assignment**, in the camera's field-angle frame: 413 = top-right
  (x>0, y>0), 439 = top-left, 465 = bottom-left, 490 = bottom-right.  This matches
  `obscuration.QUADRANT_FILTERS` in Hironao's repo (which notes the names were
  "assigned by quadrant earlier and have not been checked against the instrument").
  The handedness of the (x, y) frame on the sky (which way is East at PA = 0) is
  not pinned down here; for a 4-rotation strategy it only mirrors the pattern.
* **The MASK `NO_DATA` bit is the footprint.** Per sub-filter it flags
  (i) a cross `|x| < 5.2'` or `|y| < 5.2'`, i.e. **~10.4' wide** — wider than the
  8.18' cross used in `wide_dithering` (`obscuration.py`), presumably because the
  mask also removes the zone where light from two adjacent filters overlaps;
  (ii) the vignetted field edge beyond r ~ 46.5';
  (iii) ~2.5% of each CCD at its edges; and
  (iv) a few amplifier-sized blocks: `1_09` (413, 2 amps), `0_43` (439, 1 amp),
  `0_20` (490, 2 amps), `1_53` (490, 1 amp), a corner of `1_24` (439).
  The central CCD column (`0_12..0_15`, `1_12..1_15`) is entirely under the
  vertical arm and is not in the flats at all ("nocenter").
* **Two CCDs are missing from the 413 set:** `1_34` and `1_47` (their mirror
  images `0_34`, `0_47` are present for 465).  Whether they are dead, had no
  flat, or were dropped on purpose should be checked with Hironao.
* **The flats are in ADU, not gain-corrected.**  The steps at amplifier
  boundaries match the per-amplifier gains stored in the embedded Detector
  table exactly (e.g. `0_53` amp 0->1 step 1.265 vs gain ratio 4.41/3.53;
  `0_28` amp 3 step 0.761 vs 3.67/4.81; the header `T_GAINn` values are ~10%
  off and do *not* match).  Multiplying each amplifier by its table gain
  removes most of the CCD-to-CCD structure (scatter inside 30' drops from
  8-11% to 5-8%) and collapses the four sub-filters onto one radial curve.
* **Part of the radial decline is the optics, not the lamp.**  A flat from a
  uniform source scales with the sky area per pixel, which the HSC distortion
  reduces to 0.978 of on-axis at r = 20', 0.949 at 30', 0.903 at 40' and 0.855
  at 47' (from the embedded distortion polynomial).  A point source's flux is
  unaffected, so this Jacobian is divided out of the flat.
* **The gain- and Jacobian-corrected flat still falls with field radius**, now
  ~0.7% per arcmin from 12' to 40' (about 2/3 of the original slope), and
  steeply beyond ~43'.  Real HSC vignetting is small inside ~40', so the
  remaining gentle slope is taken to be the dome-screen illumination pattern
  and is divided out (linear fit at 12' < r < 38', extrapolated); what is left
  beyond ~43' (0.80-0.95 at 47') is treated as vignetting.  This is an
  assumption, recorded in the `RADIAL` table so it can be revisited.  Levels
  are not comparable between sub-filters (different lamp flux per narrow band),
  hence the per-filter normalisation.
* Illuminated area per sub-filter: 413 0.271 deg², 439 0.302, 465 0.301,
  490 0.291 (sum 1.165 deg² of the ~1.5 deg² HSC field).

## Hironao's `wide_dithering` (hsc-niji-survey/wide-observing-strategy)

Dither-pattern study for the same filter: HSC-SSP hex grid of pointings
(77.5' x 67.5'), each pointing observed at 4 PAs, with a small dither pattern
(recommended "4 + 2": a square of radius 23' plus an antipodal pair at 35').
It models the footprint as the real CCD rectangles (`hsc_focalplane_ccds.csv`
from obs_subaru) minus a uniform 8.18'-wide cross, 0/1 coverage, and scores the
fraction of area with n_obs <= 1 per filter plus a uniformity metric.  The
skytiling approach differs in optimising the global pointing positions on a
Fibonacci lattice for uniform coverage, and (here) in using the measured
masks/throughput instead of a model cross.

## Hironao's 4+2 pattern evaluated with the measured maps

`eval_dither_tiling -t hsc_mbq1_throughput.fits --pattern 4plus2 --plot hsc_mbq1_eval_4plus2.png`
replicates the `wide_dithering` evaluation (hex grid 77.5' x 67.5', the same
dither offsets at every pointing, each at PA 0/90/180/270, statistics over one
periodic cell) but with the per-sub-filter maps from the flats.  Fed his own
footprint model (CCD rectangles minus an 8.18' cross, no other masks) it
reproduces his numbers (phase 1 n_obs <= 1 0.48% vs 0.36%, DJS 8.30 vs 8.29;
phase 1+2 mean n_obs 6.46 vs 6.49), so the differences below come from the
footprint, not the evaluator.  Mean over the four sub-filters, 0.2' cells:

| footprint | area / quadrant | phase 1: mean n, rms/mean, n<=1, union n<=1 | phase 1+2: mean n, rms/mean, n<=1, union n<=1 |
|---|---|---|---|
| Hironao model (8.18' cross) | 0.391 deg² | 4.30, 0.25, 0.5%, 1.6% | 6.46, 0.20, 0.00%, 0.01% |
| + 10.4' cross as in the flats | 0.367 | 4.04, 0.29, 1.6%, 4.1% | 6.06, 0.22, 0.03%, 0.1% |
| + vignetted edge r > 46.5' masked | 0.308 | 3.39, 0.34, 5.8%, 10.6% | 5.09, 0.25, 0.24%, 0.8% |
| + dead amps / missing 413 CCDs (= `GOOD` maps) | 0.291 | 3.20, 0.37, 7.8%, 17.0% | 4.80, 0.28, 0.64%, 2.2% |

Per sub-filter with the `GOOD` maps, phase 1+2: 413 is the worst (mean 4.47,
rms/mean 0.32, n <= 1 1.2%) because of its two missing CCDs; 439/465/490
have mean 4.8-5.0 and n <= 1 0.3-0.6%.  Weighting by `THRU` instead of 0/1
lowers the mean depth by 2-4% and changes the uniformity negligibly, since
the maps are ~1 over the bulk.  The masked field edge is the single largest
difference from Hironao's model (it removes 15% of the area), then the wider
cross, then the dead amplifiers.  `hsc_mbq1_eval_4plus2.png` shows the
coverage maps.

### Phase 1 only (4 dithers), and tighter grids

Hironao's pointing list (`pointing/hsc_pointing_list_coadd_dr4_FCunique.csv`)
is not an all-sky tiling: it is the 802 HSC-SSP Wide pointings at
-6° < Dec < +6.4° (~1165 deg²), on rows 67.5' apart with a *fixed* RA step
of 1.299° (77.5' x cos Dec on the sky, so the spacing is only right near the
equator; at Dec 20° it would already be 73').  With only the phase-1 square
(16 exposures per pointing) and the `GOOD` maps, shrinking the whole hex grid
(dither offsets kept at their nominal values) gives:

| grid scale | spacing | exposures/deg² | mean n | rms/mean | n<=1 | n=0 | <=1 in any filter | 413 n<=1 |
|---|---|---|---|---|---|---|---|---|
| 1.00 | 77.5' x 67.5' | 11.0 | 3.20 | 0.37 | 7.8% | 1.2% | 17.0% | 11.5% |
| 0.95 | 73.6' x 64.1' | 12.2 | 3.55 | 0.35 | 5.4% | 1.4% | 12.0% | 8.8% |
| 0.90 | 69.8' x 60.8' | 13.6 | 3.95 | 0.32 | 3.7% | 0.6% | 8.0% | 5.7% |
| 0.85 | 65.9' x 57.4' | 15.2 | 4.43 | 0.31 | 1.2% | 0.05% | 4.3% | 3.5% |
| 0.80 | 62.0' x 54.0' | 17.2 | 5.00 | 0.32 | 0.7% | 0.1% | 2.4% | 1.9% |
| (4+2 at 1.00) | 77.5' x 67.5' | 16.5 | 4.80 | 0.28 | 0.6% | 0.03% | 2.2% | 1.2% |

Scaling the dither offsets along with the grid is slightly worse than keeping
them fixed.  Closing the grid removes the holes quickly but the rms/mean stays
at ~0.31-0.32: the non-uniformity is intrinsic to repeating one 16-exposure
pattern, and at equal cost (scale 0.85 vs 4+2) the 4-dither tight grid is a
little worse than 4+2 on every metric.  Improving uniformity beyond this needs
the pointing positions themselves optimised (the skytiling approach) rather
than a scaled lattice.

## Next steps

1. Confirm with Hironao: the missing 413 CCDs, the ~10.4' cross width, and
   that the illumination-gradient assumption above is reasonable.
2. Extend `optimize_tiling` to (a) take the throughput maps instead of CCD
   polygons, (b) rotate each pointing's footprint by the 4 PAs, and (c)
   optimise per-sub-filter coverage uniformity (4 maps per tile instead of 1).
3. `make_initial_tiling` for this camera: mean passes per sub-filter =
   4 x (illuminated area of that sub-filter) / tile area.
