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

## Next steps

1. Confirm with Hironao: the missing 413 CCDs, the ~10.4' cross width, and
   that the illumination-gradient assumption above is reasonable.
2. Extend `optimize_tiling` to (a) take the throughput maps instead of CCD
   polygons, (b) rotate each pointing's footprint by the 4 PAs, and (c)
   optimise per-sub-filter coverage uniformity (4 maps per tile instead of 1).
3. `make_initial_tiling` for this camera: mean passes per sub-filter =
   4 x (illuminated area of that sub-filter) / tile area.
