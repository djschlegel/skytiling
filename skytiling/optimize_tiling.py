#!/usr/bin/env python3
"""
Tiling Optimizer: All-sky survey coverage simulation optimizer.
Uses simulated annealing to minimize the RMS of coverage counts across a
catalog of random points.
"""

import argparse
import time
import sys
import math
import random
import os
import re
from typing import List, Tuple, Optional

# Limit BLAS / OpenMP threading to one thread per process BEFORE importing numpy.
# The workers do many tiny matrix ops (3x3 rotations, small dot products); with
# many worker processes (e.g. --workers 128 on a Perlmutter CPU node) each numpy
# BLAS call would otherwise spawn its own thread pool, oversubscribing the cores
# (128 processes x 128 BLAS threads) and slowing the run dramatically.  These
# must be set before numpy/BLAS import to take effect.  setdefault respects any
# value the user has already exported in the environment.
for _blas_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                  "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_blas_var, "1")

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree
import multiprocessing as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon

_shared_coverage_counts = None
_shared_randoms_xyz = None

def spherical_triangle_area(a, b, c):
    """Area of a spherical triangle with unit vertices a, b, c in steradians."""
    num = np.abs(np.dot(a, np.cross(b, c)))
    den = 1.0 + np.dot(a, b) + np.dot(b, c) + np.dot(c, a)
    return 2.0 * np.arctan2(num, den)

def is_point_in_spherical_triangle(p, a, b, c):
    """Check if unit vector p is inside the spherical triangle formed by unit vectors a, b, c."""
    # p = alpha*a + beta*b + gamma*c
    # Solve for [alpha, beta, gamma] using matrix inverse
    m = np.stack([a, b, c], axis=1)
    try:
        coeffs = np.linalg.solve(m, p)
        return np.all(coeffs > 0)
    except np.linalg.LinAlgError:
        return False

def is_point_in_ccd(p, normals, eps=1e-5):
    """Check if unit vector p is inside the spherical quad defined by inward normals.
    This version is tolerant: points whose dot product with any normal is
    less than ``-eps`` are considered outside.  It is used for coverage
    calculations where points on the boundary are acceptable.
    """
    for n in normals:
        if np.dot(p, n) < -eps:
            return False
    return True

def is_point_strictly_in_ccd(p, normals, eps=1e-8):
    """Strict interior test for overlap detection.
    Returns ``True`` only if the point lies *strictly* inside the CCD –
    i.e. the dot product with **all** inward normals is greater than ``eps``.
    This avoids treating points that lie on an edge or vertex as interior,
    which can otherwise produce false‑positive overlap reports.
    """
    # Ensure normals is a numpy array for vectorization
    normals_arr = np.asarray(normals)
    # Dot product of p with each normal (normals shape: (4,3)).
    # We require all to be > eps (strictly interior).
    return np.all(np.dot(p, normals_arr.T) > eps)


def is_point_in_ccd_robust(p, corners):
    """Check if unit vector p is inside the spherical quad defined by 4 unit corners."""
    # A convex quad can be split into two triangles: (0,1,2) and (0,2,3)
    return is_point_in_spherical_triangle(p, corners[0], corners[1], corners[2]) or \
           is_point_in_spherical_triangle(p, corners[0], corners[2], corners[3])

def arcs_intersect(v1, v2, w1, w2, eps=1e-4):
    """Check if two great-circle arcs (v1,v2) and (w1,w2) intersect strictly in their interiors."""
    n1 = np.cross(v1, v2)
    n2 = np.cross(w1, w2)
    L = np.cross(n1, n2)
    norm_L = np.linalg.norm(L)
    if norm_L < 1e-12:
        return False

    for p in [L / norm_L, -L / norm_L]:
        # p is on arc (v1, v2) if (v1 x p) . (v1 x v2) > 0 and (p x v2) . (v1 x v2) > 0
        # We use a larger epsilon to exclude endpoints.
        if (np.dot(np.cross(v1, p), n1) > eps and np.dot(np.cross(p, v2), n1) > eps and
            np.dot(np.cross(w1, p), n2) > eps and np.dot(np.cross(p, w2), n2) > eps):
            return True
    return False

def radec_to_xyz(ra, dec):
    """Convert RA/Dec (degrees) to unit XYZ vectors."""
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.stack([x, y, z], axis=-1)

def xyz_to_radec(xyz):
    """Convert unit XYZ vectors to RA/Dec (degrees)."""
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    ra = np.degrees(np.arctan2(y, x)) % 360.0
    dec = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    return ra, dec

def get_rotation_matrix(B1):
    """
    Equatorial rotation matrix that rotates the boresight [0, 0, 1] to B1.
    The resulting matrix R = [East | North | Boresight] is the Local-to-Global transformation:
    v_global = R @ v_local, where local axes are X=East, Y=North, Z=Boresight.
    This formulation avoids the singularity of the cross-product method for directions near the pole.
    """
    B1 = np.asarray(B1, dtype=np.float64)
    B1 /= np.linalg.norm(B1)

    # Handle exact poles explicitly
    if np.allclose(B1, [0, 0, 1], atol=1e-8):
        return np.eye(3)
    if np.allclose(B1, [0, 0, -1], atol=1e-8):
        # South pole: East = [1,0,0], North = [0,-1,0]
        return np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])

    # Extract RA/Dec from B1 = (cos Dec cos RA, cos Dec sin RA, sin Dec)
    x, y, z = B1
    dec = np.arcsin(np.clip(z, -1.0, 1.0))
    ra = np.arctan2(y, x)

    # East direction (local X) = -sin(ra) * i + cos(ra) * j
    east = np.array([-np.sin(ra), np.cos(ra), 0.0])
    # North direction (local Y) = -cos(ra)*sin(dec)*i - sin(ra)*sin(dec)*j + cos(dec)*k
    north = np.array([-np.cos(ra) * np.sin(dec), -np.sin(ra) * np.sin(dec), np.cos(dec)])
    # Boresight (local Z) = B1
    boresight = B1

    # Stack columns to form rotation matrix (Local -> Global)
    return np.stack([east, north, boresight], axis=1)

def parse_resume_iteration(tiles_path):
    """Detect a resume point from a tile filename of the form
    ``{something}_tiles_XXXX.fits`` where XXXX is a (zero-padded) integer.

    Returns (iteration, prefix) if the basename matches that pattern, else (None, None).
    When iteration is found the run is a continuation: the loaded tiles already
    reflect iteration XXXX, so the next iteration produced is XXXX+1.
    """
    if not tiles_path:
        return None, None
    base = os.path.basename(tiles_path)
    # Match a trailing "tiles_<digits>.fits", optionally preceded by a
    # "{something}_" prefix.
    m = re.search(r'(?:^(.+?)_)?tiles_(\d+)\.fits$', base, flags=re.IGNORECASE)
    if m:
        prefix = m.group(1) if m.group(1) else ""
        # Ensure the prefix ends with an underscore if it's not empty
        if prefix and not prefix.endswith('_'):
            prefix += '_'
        return int(m.group(2)), prefix
    return None, None


# =============================================================================
# Data Structures
# =============================================================================

class CameraGeometry:
    def __init__(self, camera_fits_path: str, verbose: bool = True):
        # Store reference rotation and boresight for later plotting
        self.R_ref = None
        self.ref_ra = None
        self.ref_dec = None
        with fits.open(camera_fits_path) as hdul:
            # Find all HDUs that contain image data
            data_hdus = []
            for hdu in hdul:
                if hdu.data is not None and hasattr(hdu.data, 'shape') and len(hdu.data.shape) == 2:
                    data_hdus.append(hdu)

            if not data_hdus:
                raise RuntimeError(f"Could not find image data in {camera_fits_path}")

            # Use CRVAL1, CRVAL2 from the first HDU's header as the reference boresight
            first_hdu = data_hdus[0]
            header = first_hdu.header
            ref_ra = header.get('CRVAL1')
            ref_dec = header.get('CRVAL2')
            if ref_ra is None or ref_dec is None:
                raise RuntimeError(f"Could not find CRVAL1 or CRVAL2 in header of {camera_fits_path}")

            center_world_ref = np.array([ref_ra, ref_dec])
            ref_boresight_xyz = radec_to_xyz(center_world_ref[0], center_world_ref[1])


            # Rotation matrix that takes [0, 0, 1] to ref_boresight_xyz
            R_ref = get_rotation_matrix(ref_boresight_xyz)
            # Store for later use and for reference boresight position (RA/Dec)
            self.R_ref = R_ref
            self.ref_ra = float(center_world_ref[0])
            self.ref_dec = float(center_world_ref[1])
            all_local_corners = []
            all_normals = []
            all_world_corners = []

            for hdu in data_hdus:
                wcs = WCS(hdu.header)
                ny, nx = hdu.data.shape

                # Corners of this CCD
                corners_pix = np.array([[0, 0], [nx, 0], [nx, ny], [0, ny]])
                corners_world = wcs.all_pix2world(corners_pix[:, 0], corners_pix[:, 1], 0)
                # Store world corners for plotting (RA/Dec in degrees)
                all_world_corners.append(np.column_stack([corners_world[0], corners_world[1]]))
                corners_xyz = radec_to_xyz(corners_world[0], corners_world[1])

                # Transform to camera local system (relative to ref_boresight)
                # local_v = R_ref^T @ global_v  =>  local_v_row = global_v_row @ R_ref
                # R_ref is the matrix that takes [0,0,1] (local) to ref_boresight (global).
                # Thus global_v = R_ref @ local_v  => local_v = R_ref.T @ global_v.
                # In matrix form: local_v_row = global_v_row @ R_ref.
                local_corners = np.dot(corners_xyz, R_ref)

                # Precompute inward local edge normals with robust orientation
                # Compute raw normals using the right-hand rule for each edge
                raw_normals = []
                for i in range(4):
                    v1 = local_corners[i]
                    v2 = local_corners[(i+1)%4]
                    n = np.cross(v1, v2)
                    raw_normals.append(n)

                # Determine the centroid of the CCD (normalized) to test orientation
                ccd_center = np.mean(local_corners, axis=0)
                norm = np.linalg.norm(ccd_center)
                if norm > 0:
                    ccd_center /= norm

                # If the first normal points outward relative to the centroid, flip all normals
                if np.dot(raw_normals[0], ccd_center) < 0:
                    raw_normals = [-n for n in raw_normals]

                normals = raw_normals

                all_local_corners.append(local_corners)
                all_normals.append(normals)

            self.ccd_corners_local = np.array(all_local_corners)
            self.normals_local = np.array(all_normals)
            self.ccd_corners_world = np.array(all_world_corners)

            # Report area and check for overlaps
            total_camera_area = 0.0
            deg2_per_sr = (180.0 / np.pi)**2
            ccd_areas = []
            if verbose:
                print("Camera Footprint Analysis:")
            for i in range(len(all_local_corners)):
                c = all_local_corners[i]
                # Split quad into two triangles: (0,1,2) and (0,2,3)
                area_sr = spherical_triangle_area(c[0], c[1], c[2]) + \
                           spherical_triangle_area(c[0], c[2], c[3])
                area_deg2 = area_sr * deg2_per_sr
                total_camera_area += area_deg2
                ccd_areas.append(area_deg2)
                if verbose:
                    print(f"  CCD {i:2}: {area_deg2:8.4f} sq deg")
            if verbose:
                print(f"  Total Camera Area: {total_camera_area:8.4f} sq deg\n")
            # On-sky area of each CCD and of the whole footprint [deg^2].
            # Used by make_initial_tiling to convert a mean number of passes
            # into a Fibonacci tile area.
            self.ccd_areas_deg2 = np.array(ccd_areas)
            self.total_area_deg2 = float(total_camera_area)

            # Overlap check (always performed; report overlaps without stopping execution)
            overlap_issues = []
            for i in range(len(all_local_corners)):
                for j in range(i + 1, len(all_local_corners)):
                    # Check for interior overlaps.
                    # Note: shared boundaries/corners are NOT overlaps.
                    # We check if any corner of one CCD is strictly inside the other.
                    overlap_found = False

                    # 1. Check corners of i in j
                    for p in all_local_corners[i]:
                        if is_point_strictly_in_ccd(p, all_normals[j], eps=1e-8):
                            overlap_found = True
                            reason = f"corner of {i} in {j}"
                            break

                    if not overlap_found:
                        # 2. Check corners of j in i
                        for p in all_local_corners[j]:
                            if is_point_strictly_in_ccd(p, all_normals[i], eps=1e-8):
                                overlap_found = True
                                reason = f"corner of {j} in {i}"
                                break

                    if not overlap_found:
                        # 3. Check edge intersections
                        for ei in range(4):
                            v1, v2 = all_local_corners[i][ei], all_local_corners[i][(ei+1)%4]
                            for ej in range(4):
                                w1, w2 = all_local_corners[j][ej], all_local_corners[j][(ej+1)%4]
                                if arcs_intersect(v1, v2, w1, w2):
                                    overlap_found = True
                                    reason = "edge intersection"
                                    break
                            if overlap_found:
                                break

                    if overlap_found:
                        overlap_issues.append((i, j, reason))

            if overlap_issues:
                print("WARNING: Overlapping CCDs detected:")
                for i, j, reason in overlap_issues:
                    print(f"  CCD {i} and CCD {j}: {reason}")

            # max_radius = max angle from [0,0,1] to any corner of any CCD
            # The z-component of local_corners is dot(local_corner, [0,0,1])
            all_z = self.ccd_corners_local[:, :, 2]
            self.max_radius = np.max(np.arccos(np.clip(all_z, -1, 1)))
            if verbose:
                print(f"Camera max_radius: {np.degrees(self.max_radius):.4f} degrees")
            self.R_ref = R_ref

class RandomCatalog:
    def __init__(self, fits_path: Optional[str], num_randoms: int):
        if fits_path:
            with fits.open(fits_path) as hdul:
                data = hdul[1].data
                ra = data['RA'] if 'RA' in data.names else data.field(0)
                dec = data['DEC'] if 'DEC' in data.names else data.field(1)
                xyz = radec_to_xyz(ra, dec).astype(np.float64)
        else:
            print(f"Generating {num_randoms} random points...")
            phi = np.random.uniform(0, 2 * np.pi, num_randoms)
            costheta = np.random.uniform(-1, 1, num_randoms)
            sintheta = np.sqrt(1 - costheta**2)
            xyz = np.stack([
                sintheta * np.cos(phi),
                sintheta * np.sin(phi),
                costheta
            ], axis=-1).astype(np.float64)

            # Write generated randoms to a FITS file for reproducibility
            ra, dec = xyz_to_radec(xyz)
            col_ra = fits.Column(name='RA', format='D', array=ra)
            col_dec = fits.Column(name='DEC', format='D', array=dec)
            hdu = fits.BinTableHDU.from_columns([col_ra, col_dec])
            filename = f"generated_randoms_{num_randoms}.fits"
            hdu.writeto(filename, overwrite=True)
            print(f"Generated randoms saved to {filename}")

        self.n_rand = len(xyz)

        # Back the randoms coordinates with a shared RawArray so worker
        # processes inherit it via fork (no per-task pickling of the full
        # array).  Shape is (n_rand, 3) flattened in C order.
        xyz = np.ascontiguousarray(xyz, dtype=np.float64)
        self._raw_xyz = mp.RawArray('d', self.n_rand * 3)
        self.xyz = np.frombuffer(self._raw_xyz, dtype=np.float64).reshape(self.n_rand, 3)
        self.xyz[:] = xyz

        self._raw_counts = mp.RawArray('i', self.n_rand)
        self.coverage_counts = np.frombuffer(self._raw_counts, dtype=np.int32)

class TileState:
    def __init__(self, initial_xyz: np.ndarray, fixed: bool = False):
        self.center_xyz = initial_xyz.copy()
        self.initial_xyz = initial_xyz.copy()
        self.cand_idx = np.array([], dtype=np.int32)
        self.fixed = fixed

# =============================================================================
# Worker Logic
# =============================================================================

def _inside_any_ccd(pts_local, camera_normals):
    """Vectorized inside-any-CCD test.

    pts_local: (M, 3) candidate points already rotated into the tile frame.
    camera_normals: (n_ccd, 4, 3) inward edge normals per CCD.
    Returns a boolean array (M,) that is True where a point is inside ANY CCD
    (i.e. on the positive side of all 4 edge normals of at least one CCD).
    """
    # dots: (M, n_ccd, 4) = point . each edge normal of each CCD.
    dots = np.einsum('md,ced->mce', pts_local, camera_normals)
    # inside a given CCD => all 4 edge dots >= 0; inside any CCD => any over CCDs.
    return np.any(np.all(dots >= 0.0, axis=2), axis=1)


def worker_optimize_subset(subset_tiles_data,
                           n_rand,
                           sum_c,
                           camera_normals,
                           sigma,
                           temp,
                           max_drift,
                           seed,
                           moves_per_tile=1):
    """
    Worker function to process a subset of tiles in a layer.

    For each tile up to ``moves_per_tile`` candidate moves are proposed from the
    tile's current position.  The inside/outside state at the current position
    (``inside_old``) is computed once and reused for every proposal, and the
    per-candidate inside test is vectorized across all candidates.  The proposal
    with the lowest energy change (best descent) is selected, then run through
    the Metropolis acceptance test.

    Returns: list of (tile_id, new_xyz, accepted, improve_worsen, delta_sum_c, delta_sum_c2, move_dist)
    """
    rng = np.random.default_rng(seed)

    # Access the shared arrays from the global scope (inherited via fork).
    # Neither the randoms coordinates nor the coverage counts are pickled per
    # task; workers attach to the parent's RawArrays directly.
    global _shared_coverage_counts, _shared_randoms_xyz
    coverage_counts = np.frombuffer(_shared_coverage_counts, dtype=np.int32)
    randoms_xyz = np.frombuffer(_shared_randoms_xyz, dtype=np.float64).reshape(n_rand, 3)

    # Running global sum of coverage counts (Sigma c).  We minimize the true
    # variance N*Var = Sigma c^2 - (Sigma c)^2 / N, which is monotonic in RMS,
    # rather than Sigma c^2 alone (which ignores mean changes and can accept
    # moves that increase RMS).  Tiles in a layer have disjoint candidate sets,
    # so updating this running value as we accept moves within the subset keeps
    # it exact for this worker's sequential tiles.
    running_sum_c = float(sum_c)

    results = []

    for tile_id, current_xyz, initial_xyz, cand_idx in subset_tiles_data:
        if len(cand_idx) == 0:
            results.append((tile_id, current_xyz, False, 0, 0, 0.0, 0.0))
            continue

        # Candidate point coordinates and their current coverage counts.
        pts = randoms_xyz[cand_idx]                       # (C, 3)
        cur_counts = coverage_counts[cand_idx].astype(np.int64)  # (C,)

        # --- inside_old: computed ONCE per tile (independent of the proposal) ---
        R_old = get_rotation_matrix(current_xyz)
        inside_old = _inside_any_ccd(pts @ R_old, camera_normals)  # (C,) bool

        # Local frame axes for proposing tangential offsets (constant per tile).
        z_axis = current_xyz
        x_axis = np.array([1, 0, 0]) if abs(z_axis[0]) < 0.9 else np.array([0, 1, 0])
        x_axis = np.cross(x_axis, z_axis)
        x_axis /= np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)

        # --- Evaluate up to moves_per_tile proposals; keep the best descent ---
        best = None  # (delta_energy, proposed_xyz, delta_sum_c, delta_sum_c2, changed_mask, inside_new)
        for _ in range(max(1, int(moves_per_tile))):
            dx = rng.normal(0, sigma)
            dy = rng.normal(0, sigma)
            proposed_xyz = current_xyz + dx * x_axis + dy * y_axis
            proposed_xyz /= np.linalg.norm(proposed_xyz)

            dist_from_init = np.arccos(np.clip(np.dot(proposed_xyz, initial_xyz), -1, 1))
            if dist_from_init > max_drift:
                proposed_xyz = initial_xyz + (proposed_xyz - initial_xyz) * (max_drift / dist_from_init)
                proposed_xyz /= np.linalg.norm(proposed_xyz)

            R_new = get_rotation_matrix(proposed_xyz)
            inside_new = _inside_any_ccd(pts @ R_new, camera_normals)  # (C,) bool

            changed = inside_old != inside_new
            if not np.any(changed):
                continue

            # +1 where a point entered coverage, -1 where it left.
            delta = np.where(inside_new[changed], 1, -1).astype(np.int64)
            old_c = cur_counts[changed]
            new_c = old_c + delta
            d_sum_c = int(delta.sum())
            d_sum_c2 = float(np.sum(new_c**2 - old_c**2))

            # Energy change (N*Var); minimizing this minimizes RMS.
            d_energy = d_sum_c2 - (2.0 * running_sum_c * d_sum_c + d_sum_c**2) / n_rand

            if best is None or d_energy < best[0]:
                best = (d_energy, proposed_xyz, d_sum_c, d_sum_c2, changed, inside_new)

        if best is None:
            # No proposal changed any candidate's coverage state.
            results.append((tile_id, current_xyz, False, 0, 0, 0.0, 0.0))
            continue

        delta_energy, proposed_xyz, delta_sum_c, delta_sum_c2, changed, inside_new = best

        # --- Metropolis acceptance on the best proposal ---
        accepted = False
        improve_worsen = 0
        if delta_energy < 0:
            accepted = True
            improve_worsen = 1
        elif temp > 0 and rng.random() < np.exp(-delta_energy / (temp * 1e6)):
            accepted = True
            improve_worsen = -1

        if accepted:
            # Apply +/-1 to the coverage counts of the candidates that changed.
            changed_idx = cand_idx[changed]
            coverage_counts[changed_idx] += np.where(inside_new[changed], 1, -1).astype(np.int32)
            running_sum_c += delta_sum_c

            # Distance between old and new positions (in degrees).
            move_dist = np.degrees(np.arccos(np.clip(np.dot(current_xyz, proposed_xyz), -1, 1)))
            results.append((tile_id, proposed_xyz, True, improve_worsen, delta_sum_c, delta_sum_c2, move_dist))
        else:
            results.append((tile_id, current_xyz, False, 0, 0, 0.0, 0.0))

    return results

# =============================================================================
# Main Optimizer
# =============================================================================

class TilingOptimizer:
    def __init__(self, camera_fits, tiles_fits, randoms_fits,
                 num_randoms=10_000_000, max_drift=0.3,
                 delta=0.1, shrink=None, temp=0.1,
                 iters=1000, workers=None, moves_per_tile=1,
                 output_prefix="output_",
                 write_random_coverage_hdu=False,
                 plot=False, ra_center=None, dec_center=None, diameter=None,
                 fixed_ra_min=None, fixed_ra_max=None,
                 fixed_dec_min=None, fixed_dec_max=None):

        _t_init_start = time.time()

        _t0 = time.time()
        self.camera = CameraGeometry(camera_fits)
        print(f"[timing] Camera geometry loaded in {time.time() - _t0:.3f} s")
        # Propagate camera reference boresight for diagnostic plots
        self.ref_ra = self.camera.ref_ra
        self.ref_dec = self.camera.ref_dec
        _t0 = time.time()
        self.randoms = RandomCatalog(randoms_fits, num_randoms)
        print(f"[timing] Random catalog prepared in {time.time() - _t0:.3f} s")
        self.write_random_coverage_hdu = write_random_coverage_hdu
        self.plot = plot
        self.ra_center = ra_center
        self.dec_center = dec_center
        self.diameter = diameter
        # Store optional bounding box limits for fixing tiles
        self.fixed_ra_min = fixed_ra_min
        self.fixed_ra_max = fixed_ra_max
        self.fixed_dec_min = fixed_dec_min
        self.fixed_dec_max = fixed_dec_max

        # Open the input tiles FITS file and store its column definitions and original data for later reuse.
        _t0 = time.time()
        with fits.open(tiles_fits) as hdul:
            # Preserve the original column definitions (names, formats, etc.)
            self._orig_coldefs = hdul[1].columns
            # Preserve the original data arrays for all columns
            self._orig_data = {col.name: np.array(hdul[1].data[col.name]) for col in self._orig_coldefs}

            data = hdul[1].data
            # Map UPPERCASE column name -> actual column name for case-insensitive
            # access (e.g. the all-sky files use lowercase 'ra'/'dec'/'tileid').
            name_by_upper = {n.upper(): n for n in data.names}
            ra = data[name_by_upper['RA']] if 'RA' in name_by_upper else data.field(0)
            dec = data[name_by_upper['DEC']] if 'DEC' in name_by_upper else data.field(1)
            # Optional FIXED column: 1 = fixed, 0 = free (default 0)
            if 'FIXED' in name_by_upper:
                fixed_col = data[name_by_upper['FIXED']]
            else:
                fixed_col = np.zeros(len(ra), dtype=int)
            initial_xyz = radec_to_xyz(ra, dec)

            # Determine fixed flags based on column and bounding boxes
            tile_ra, tile_dec = xyz_to_radec(initial_xyz)
            fixed_flags = (fixed_col == 1)
            # Apply external bounding boxes (if provided)
            if self.fixed_ra_min is not None:
                fixed_flags |= (tile_ra < self.fixed_ra_min)
            if self.fixed_ra_max is not None:
                fixed_flags |= (tile_ra > self.fixed_ra_max)
            if self.fixed_dec_min is not None:
                fixed_flags |= (tile_dec < self.fixed_dec_min)
            if self.fixed_dec_max is not None:
                fixed_flags |= (tile_dec > self.fixed_dec_max)

        self.tiles = [TileState(xyz, fixed=bool(fixed_flags[i])) for i, xyz in enumerate(initial_xyz)]
        # Diagnostic information about the tile set
        num_tiles = len(self.tiles)
        print(f"Total tiles loaded: {num_tiles}")
        # tile_ra and tile_dec were computed just above from initial_xyz
        print(f"RA range of tiles: {tile_ra.min():.3f} - {tile_ra.max():.3f} deg")
        print(f"Dec range of tiles: {tile_dec.min():.3f} - {tile_dec.max():.3f} deg")
        # Report how many tiles ended up fixed (for debugging)
        fixed_count = sum(t.fixed for t in self.tiles)
        print(f"Fixed tiles after applying bounds: {fixed_count}")
        print(f"[timing] Tile catalog loaded in {time.time() - _t0:.3f} s")
        self.n_tiles = len(self.tiles)

        self.max_drift = np.radians(max_drift)
        self.sigma0 = np.radians(delta)
        if shrink is not None:
            self.shrink = shrink
        elif iters > 0:
            self.shrink = (iters - 1) / iters
        else:
            self.shrink = 1.0
        self.temp = temp
        self.iters = iters
        self.workers = workers or mp.cpu_count()
        self.moves_per_tile = max(1, int(moves_per_tile))
        # Prefix applied to every output file (FITS checkpoints and PNG plots).
        self.output_prefix = output_prefix

        # If the input tile file looks like a previous checkpoint
        # ({something}_tiles_XXXX.fits), resume after iteration XXXX so the next
        # files written are XXXX+1, etc.
        res_iter, res_prefix = parse_resume_iteration(tiles_fits)
        if res_iter is not None:
            self.start_iter = res_iter
            # If a prefix was detected in the filename, it takes precedence.
            if res_prefix:
                self.output_prefix = res_prefix
            print(f"Resuming from iteration {res_iter} (input '{os.path.basename(tiles_fits)}'); "
                  f"next output iteration is {res_iter + 1}.")
            if res_prefix:
                print(f"Using resumed output prefix: {self.output_prefix}")
        else:
            self.start_iter = 0

        _t0 = time.time()
        self._precompute_candidates()
        print(f"[timing] Candidate precompute in {time.time() - _t0:.3f} s")

        print("Computing initial coverage...")
        _t0 = time.time()
        self.randoms.coverage_counts.fill(0)
        for i, tile in enumerate(self.tiles):
            if i % 500 == 0:
                print(f"  Processing tile {i}/{self.n_tiles}...")

            R = get_rotation_matrix(tile.center_xyz)
            n_local_all = self.camera.normals_local
            n_global = n_local_all @ R.T

            if len(tile.cand_idx) == 0:
                continue
            pts = self.randoms.xyz[tile.cand_idx]  # (N_cand, 3)

            # Check if point is inside ANY CCD
            # A point is inside a CCD if it's on the positive side of all 4 normals
            inside_any = np.zeros(len(pts), dtype=bool)
            for ccd_normals in n_global:
                inside_this_ccd = np.ones(len(pts), dtype=bool)
                for n in ccd_normals:
                    inside_this_ccd &= (pts @ n >= 0)
                    if not np.any(inside_this_ccd):
                        break
                inside_any |= inside_this_ccd
                if np.all(inside_any):
                    break

            self.randoms.coverage_counts[tile.cand_idx[inside_any]] += 1

        print(f"[timing] Initial coverage computed in {time.time() - _t0:.3f} s")
        counts = self.randoms.coverage_counts
        self.sum_c = np.sum(counts).astype(np.float64)
        self.sum_c2 = np.sum(counts.astype(np.float64)**2)
        self.mean_coverage = self.sum_c / self.randoms.n_rand
        print(f"Mean coverage = {self.mean_coverage:.4f}")

        # Plot camera footprint if requested.  This is the only place it is
        # produced ({prefix}camera.png); it depends only on the boresight, not
        # on iteration, so it has no numeric suffix.  The per-iteration coverage
        # and layers plots are written by _save_iteration_artifacts (starting at
        # iteration 0), so we do not duplicate them here.
        if self.plot:
            self._plot_camera_footprint()
        _t0 = time.time()
        self.layers = self._build_layers()
        print(f"[timing] Layer partitioning in {time.time() - _t0:.3f} s")
        print(f"Tiling partitioned into {len(self.layers)} layers.")
        # Report initial RMS (before any optimization steps)
        initial_rms = self.get_rms()
        print(f"Initial RMS = {initial_rms:.10f}")

        print(f"[timing] Total initialization time: {time.time() - _t_init_start:.3f} s")

    def _precompute_candidates(self):
        print("Precomputing candidate indices...")
        tree = cKDTree(self.randoms.xyz)
        search_radius = self.camera.max_radius + self.max_drift
        for tile in self.tiles:
            tile.cand_idx = np.array(tree.query_ball_point(tile.center_xyz, search_radius), dtype=np.int32)

    def _build_layers(self) -> List[List[int]]:
        """Partition tiles into layers such that no two tiles within angular
        distance ``sep`` share a layer (so their candidate-random sets are
        disjoint and workers never write the same coverage index).

        This is a graph-coloring problem.  We use a cKDTree to find, for each
        tile, only the handful of tiles within ``sep`` (its conflict
        neighbours), then greedily color in first-fit-decreasing order
        (descending candidate count, same heuristic as before).  This is
        near-linear in the number of tiles rather than the previous O(N^2)
        all-pairs scan, which slowed dramatically for large all-sky grids.
        """
        sep = 2 * (self.camera.max_radius + self.max_drift)
        print(f"Partitioning {self.n_tiles} tiles into layers...")

        if self.n_tiles == 0:
            return []

        # Tile centers as unit vectors; build a spatial index over them.
        pts = np.array([t.initial_xyz for t in self.tiles])
        tree = cKDTree(pts)

        # Convert the angular separation threshold to a Euclidean chord radius
        # between unit vectors: chord = 2*sin(theta/2) for theta in [0, pi].
        # (For sep >= pi every pair is within range, so use the max chord 2.0.)
        if sep >= np.pi:
            chord_r = 2.0
        else:
            chord_r = 2.0 * np.sin(sep / 2.0)

        # Conflict neighbours within sep for every tile, computed in parallel.
        # neighbours[i] is a list of tile indices (including i itself).
        # Note: this uses chord <= chord_r (i.e. angular <= sep), vs the old
        # strict dist < sep; the difference is only the measure-zero boundary
        # and errs toward slightly larger separation, which is always safe.
        neighbours = tree.query_ball_point(pts, chord_r, workers=self.workers)

        # First-fit-decreasing graph coloring: same processing order as before.
        tile_indices = sorted(range(self.n_tiles),
                              key=lambda i: len(self.tiles[i].cand_idx), reverse=True)

        colors = np.full(self.n_tiles, -1, dtype=np.int64)  # -1 = unassigned
        n_layers = 0
        for count, idx in enumerate(tile_indices):
            if count % 20000 == 0:
                print(f"  Processing tile {count}/{self.n_tiles}...")
            # Colors already taken by this tile's (already-assigned) neighbours.
            forbidden = set()
            for nb in neighbours[idx]:
                c = colors[nb]
                if c >= 0:
                    forbidden.add(c)
            # Smallest layer index not used by any conflicting neighbour.
            color = 0
            while color in forbidden:
                color += 1
            colors[idx] = color
            if color + 1 > n_layers:
                n_layers = color + 1

        # Materialize layer membership lists, preserving processing order.
        layers = [[] for _ in range(n_layers)]
        for idx in tile_indices:
            layers[colors[idx]].append(idx)
        return layers

    def get_rms(self):
        var = (self.sum_c2 / self.randoms.n_rand) - (self.mean_coverage**2)
        return np.sqrt(max(0, var))

    def optimize(self):
        # Ensure workers have access to the shared memory arrays without pickling them.
        # We use module-level globals (inherited via fork) for the workers to attach to.
        global _shared_coverage_counts, _shared_randoms_xyz
        _shared_coverage_counts = self.randoms._raw_counts
        _shared_randoms_xyz = self.randoms._raw_xyz

        # Save a snapshot of the starting state.  For a fresh run this is
        # iteration 0.  For a resumed run the start snapshot already exists as
        # the input file ({prefix}tiles_XXXX.fits), so we skip re-writing it to
        # avoid clobbering the input; only genuinely new iterations are written.
        if self.start_iter == 0:
            self._save_iteration_artifacts(0)
        else:
            print(f"Resumed run: start snapshot (iteration {self.start_iter}) "
                  f"already exists; not re-writing it.")

        if self.iters == 0:
            print("iters=0: Optimization skipped. Initial artifacts saved.")
            return

        pool = mp.Pool(self.workers)
        start_time = time.time()

        # Precompute, per layer, the list of FREE (non-fixed) tile ids.  Fixed
        # tiles never move, so excluding them up front means each work chunk
        # carries a balanced number of tiles that actually do work — otherwise
        # spatially-clustered fixed tiles (e.g. from --fixed-dec-min/max bands)
        # would leave some chunks nearly empty and others full, wasting workers.
        free_by_layer = [
            [tid for tid in layer if not getattr(self.tiles[tid], 'fixed', False)]
            for layer in self.layers
        ]

        try:
            for i in range(self.iters):
                _t_iter_start = time.time()
                sigma = self.sigma0 * (self.shrink**i)
                acc_total, acc_improve, acc_worsen, rej = 0, 0, 0, 0
                accepted_dists = []

                for layer_idx, free_ids in enumerate(free_by_layer):
                    if not free_ids:
                        continue
                    # Chunk the free tiles evenly across workers so every chunk
                    # has roughly the same amount of real work.
                    subset_size = math.ceil(len(free_ids) / self.workers)
                    futures = []
                    for j in range(0, len(free_ids), subset_size):
                        subset_ids = free_ids[j : j + subset_size]
                        subset_data = [
                            (tid, self.tiles[tid].center_xyz, self.tiles[tid].initial_xyz, self.tiles[tid].cand_idx)
                            for tid in subset_ids
                        ]
                        # Generate a unique seed for each worker call
                        seed = random.randint(0, 2**31)
                        futures.append(pool.apply_async(
                            worker_optimize_subset,
                            args=(subset_data, self.randoms.n_rand, self.sum_c,
                                  self.camera.normals_local, sigma, self.temp, self.max_drift, seed,
                                  self.moves_per_tile)
                        ))

                    for f in futures:
                        res = f.get()
                        for tid, new_xyz, accepted, imp_wors, delta_c, delta_c2, m_dist in res:
                            if accepted:
                                self.tiles[tid].center_xyz = new_xyz
                                self.sum_c += delta_c
                                self.sum_c2 += delta_c2
                                self.mean_coverage = self.sum_c / self.randoms.n_rand
                                acc_total += 1
                                if imp_wors == 1: acc_improve += 1
                                else: acc_worsen += 1
                                # track accepted move distance
                                accepted_dists.append(m_dist)
                            else:
                                rej += 1
                    la = len(free_ids) # unused but keeping current context
                    la_idx = layer_idx + 1
                    la_total = len(self.layers)
                    print(f"  Layer {la_idx}/{la_total} complete", end='\r')
                    sys.stdout.flush()

                iter_time = time.time() - _t_iter_start
                rms = self.get_rms()
                rate = acc_total / (acc_total + rej) if (acc_total + rej) > 0 else 0

                if accepted_dists:
                    min_dist = np.min(accepted_dists)
                    max_dist = np.max(accepted_dists)
                    mean_dist = np.mean(accepted_dists)
                    dist_str = f" | MoveDist = [{min_dist:.4f}, {mean_dist:.4f}, {max_dist:.4f}]°"
                else:
                    dist_str = " | MoveDist = N/A"

                # Absolute iteration number, offset by the resume point.
                iter_no = self.start_iter + i + 1
                print(f"Iter {iter_no:4} | RMS = {rms:.10f} | σ = {np.degrees(sigma):.4f}° | "
                      f"Accepted = {acc_total:4} (imp={acc_improve}, wor={acc_worsen}) | "
                      f"Rejected = {rej:4} | Rate = {rate:.2%}{dist_str} | Time = {iter_time:.3f} s")

                # Save artifacts for this iteration (start_iter + i + 1 so a
                # resumed run continues the numbering from the input file).
                self._save_iteration_artifacts(iter_no)

                # Optional checkpoint at regular intervals removed per user request; checkpointing is handled by per‑iteration artifact saving.

        finally:
            pool.close()
            pool.join()

    def save_checkpoint(self, path, iteration):
        # Optimized RA/DEC from the current tile positions.
        ra_opt, dec_opt = xyz_to_radec(np.array([t.center_xyz for t in self.tiles]))

        # Locate the original RA/DEC columns case-insensitively (the all-sky
        # files use lowercase 'ra'/'dec').  We update these in place so the
        # output never ends up with a duplicate set of RA/DEC columns.
        upper_to_name = {c.name.upper(): c.name for c in self._orig_coldefs}
        ra_name = upper_to_name.get('RA')
        dec_name = upper_to_name.get('DEC')

        # Build the output RA/DEC. For fixed tiles, restore the exact original
        # values so a fixed tile is always identical to its input coordinates.
        fixed_mask = np.array([getattr(t, 'fixed', False) for t in self.tiles])
        out_ra = ra_opt.copy()
        out_dec = dec_opt.copy()
        if ra_name is not None:
            out_ra[fixed_mask] = self._orig_data[ra_name][fixed_mask]
        if dec_name is not None:
            out_dec[fixed_mask] = self._orig_data[dec_name][fixed_mask]

        # Primary HDU (empty) required by astropy
        primary_hdu = fits.PrimaryHDU()

        # Construct columns for the tile table, updating RA/DEC in place.
        cols = []
        for col in self._orig_coldefs:
            name = col.name
            fmt = col.format

            if name == ra_name:
                data = out_ra
            elif name == dec_name:
                data = out_dec
            else:
                data = self._orig_data[name]

            cols.append(fits.Column(name=name, format=fmt, array=data))

        # If RA or DEC were not present in the original file, add a single
        # column each (fallback only — avoids creating duplicates above).
        if ra_name is None:
            cols.append(fits.Column(name='RA', format='D', array=out_ra))
        if dec_name is None:
            cols.append(fits.Column(name='DEC', format='D', array=out_dec))

        hdu_tiles = fits.BinTableHDU.from_columns(cols)
        hdu_tiles.header['ITER'] = iteration
        hdu_tiles.header['DATE'] = time.strftime('%Y-%m-%dT%H:%M:%S')

        # Assemble HDU list starting with primary
        hdulist = fits.HDUList([primary_hdu, hdu_tiles])
        if self.write_random_coverage_hdu:
            # Second HDU with random point coverage counts, row‑matched to randoms file
            counts = self.randoms.coverage_counts.astype(np.int32)
            col_counts = fits.Column(name='coverage', format='J', array=counts)
            hdu_counts = fits.BinTableHDU.from_columns([col_counts])
            hdu_counts.header['COMMENT'] = 'Coverage counts for each random point'
            hdulist.append(hdu_counts)

        hdulist.writeto(path, overwrite=True)
        print(f"Checkpoint saved to {path} at iteration {iteration}")

    def _plot_camera_footprint(self):
        """Plot the CCD outlines of the camera centered at its boresight."""
        print(f"Generating camera footprint plot centered at RA={self.ref_ra:.4f}, Dec={self.ref_dec:.4f}...")

        # Use a simple gnomic projection centered at the boresight
        # We transform everything to local coordinates, then just plot (x, y)
        plt.figure(figsize=(8, 8))
        ax = plt.gca()

        # The world coordinates (RA/Dec) are stored in self.camera.ccd_corners_world
        for i, world_corners in enumerate(self.camera.ccd_corners_world):
            # world_corners shape: (4, 2) with columns [RA, Dec] in degrees
            # Compute offsets from the reference boresight, scaling RA by cos(dec)
            delta_ra = (world_corners[:, 0] - self.ref_ra) * np.cos(np.radians(self.ref_dec))
            delta_dec = world_corners[:, 1] - self.ref_dec
            pts = np.column_stack([delta_ra, delta_dec])
            # Close the polygon
            poly_pts = np.vstack([pts, pts[0]])
            ax.plot(poly_pts[:, 0], poly_pts[:, 1], label=f"CCD {i}")

            # Label each CCD in the center
            centroid = np.mean(pts, axis=0)
            ax.text(centroid[0], centroid[1], str(i),
                    ha='center', va='center', fontsize=8, color='black')

        ax.set_aspect('equal')
        ax.set_xlabel("Local X (East) [unit vectors]")
        ax.set_ylabel("Local Y (North) [unit vectors]")
        ax.set_title(f"Camera Footprint at Boresight\n(RA={self.ref_ra:.4f}, Dec={self.ref_dec:.4f})")
        # Labels are now added inside the loop for all CCDs
        if len(self.camera.ccd_corners_local) < 20:
            ax.legend()

        plt.grid(True, alpha=0.3)
        out_path = f"{self.output_prefix}camera.png"
        plt.savefig(out_path)
        plt.close()
        print(f"Saved {out_path}")
        # End of camera footprint plot; other diagnostics are handled separately
        return

    def _plot_allsky_layer0(self, out_path=None):
        """Plot all tiles in the first layer on an Aitoff all‑sky projection.
        The projection is centered at RA=120° (RA runs right‑to‑left).

        ``out_path`` is the file written directly (no intermediate rename).
        Defaults to ``{output_prefix}layers.png`` when not given.
        """
        if out_path is None:
            out_path = f"{self.output_prefix}layers.png"
        print("Generating all‑sky Aitoff plot of tiles in layer 0...")
        ra_center_deg = 120.0

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111, projection='aitoff')

        # Tiles that belong to the first layer.
        layer0_indices = self.layers[0]
        ra, dec = xyz_to_radec(np.array([self.tiles[i].center_xyz for i in layer0_indices]))

        # Shift RA to put ra_center_deg at the centre and wrap into [-180, 180]°.
        ra_shift = ra - ra_center_deg
        ra_shift = ((ra_shift + 180) % 360) - 180

        # Convert to radians for the Aitoff projection.
        ra_rad = np.radians(ra_shift)
        dec_rad = np.radians(dec)

        # Scatter the tile centers.
        ax.scatter(ra_rad, dec_rad, s=1, color='blue', alpha=0.5)

        # Custom RA formatter to convert relative radians back to absolute RA degrees.
        def format_ra(x, pos):
            ra_val = np.degrees(x) + ra_center_deg
            return f"{int(ra_val % 360)}°"

        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_ra))

        # RA increases to the left in astronomy.
        # Since invert_xaxis() is not supported for Aitoff, we can't easily flip the axes
        # but the coordinates are already in the correct relative order for the projection.
        # To strictly adhere to RA increasing to the left, the projection normally handles this
        # if the coordinates are correct. We'll remove the forbidden call.
        pass

        ax.grid(True)
        ax.set_title(f"All‑sky tiling (Layer 0) – Aitoff projection (RA centre {ra_center_deg}°)")
        plt.savefig(out_path)
        plt.close()
        print(f"Saved {out_path}")

    def plot_diagnostics(self, out_path=None):
        """Generate diagnostic plots for the tiling configuration as a coverage heatmap.

        ``out_path`` is the file written directly (no intermediate rename).
        Defaults to ``{output_prefix}coverage.png`` when not given.
        """
        if out_path is None:
            out_path = f"{self.output_prefix}coverage.png"
        if self.ra_center is None or self.dec_center is None:
            print("Skipping diagnostic plot: ra_center or dec_center not provided.")
            return

        if self.diameter is None:
            self.diameter = 5.0  # Default zoom diameter

        print(f"Generating tiling diagnostics coverage map centered at RA={self.ra_center:.4f}, Dec={self.dec_center:.4f}...")

        # Grid resolution for the heatmap
        res = 300
        # Create a grid of local offsets (degrees)
        x = np.linspace(-self.diameter / 2, self.diameter / 2, res)
        y = np.linspace(-self.diameter / 2, self.diameter / 2, res)
        X, Y = np.meshgrid(x, y)

        # Convert grid points to RA/Dec and then to XYZ unit vectors
        # Scale X by cos(dec_center) to get actual RA offsets
        cos_dec_center = np.cos(np.radians(self.dec_center))
        grid_ra = self.ra_center + X / cos_dec_center
        grid_dec = self.dec_center + Y

        # Flatten for easier vectorization
        grid_xyz = radec_to_xyz(grid_ra.flatten(), grid_dec.flatten()) # (res*res, 3)
        coverage = np.zeros(res * res, dtype=np.int32)

        # To optimize, only check tiles that are close enough to the plot center
        center_xyz = radec_to_xyz(self.ra_center, self.dec_center)
        half_side_rad = np.radians(self.diameter / 2)
        half_diag_rad = half_side_rad * np.sqrt(2)
        max_reach = half_diag_rad + self.camera.max_radius

        relevant_tiles = []
        for tile in self.tiles:
            dist = np.arccos(np.clip(np.dot(tile.center_xyz, center_xyz), -1, 1))
            if dist < max_reach:
                relevant_tiles.append(tile)

        # Calculate coverage for each grid point
        for tile in relevant_tiles:
            R = get_rotation_matrix(tile.center_xyz)
            # Transform all grid points to this tile's local system: p_local = R.T @ p_global
            # grid_xyz: (N, 3), R: (3, 3) -> p_local: (N, 3)
            p_local = grid_xyz @ R

            # For each CCD, check if point is inside (all normals positive)
            inside_any_ccd = np.zeros(len(grid_xyz), dtype=bool)
            for ccd_normals in self.camera.normals_local:
                # ccd_normals: (4, 3), p_local: (N, 3) -> (N, 4)
                dots_ccd = p_local @ ccd_normals.T
                inside_this_ccd = np.all(dots_ccd >= 0, axis=1)
                inside_any_ccd |= inside_this_ccd

            coverage += inside_any_ccd.astype(np.int32)

        # Reshape and plot
        coverage_map = coverage.reshape((res, res))

        plt.figure(figsize=(8, 8))
        im = plt.imshow(coverage_map, extent=[-self.diameter/2, self.diameter/2, -self.diameter/2, self.diameter/2],
                        origin='lower', cmap='viridis', aspect='equal', interpolation='nearest')

        plt.colorbar(im, label='Coverage Count')
        plt.xlabel("Delta RA * cos(Dec) [deg]")
        plt.ylabel("Delta Dec [deg]")
        plt.title(f"Tiling Coverage Heatmap (Zoom)\nCenter: RA={self.ra_center:.4f}, Dec={self.dec_center:.4f}, Dia={self.diameter:.2f}")
        plt.grid(False)
        plt.savefig(out_path)
        plt.close()
        print(f"Saved {out_path}")

    def _save_iteration_artifacts(self, iteration: int):
        """Save FITS checkpoint and diagnostic PNGs for a given iteration.

        * iteration 0 → snapshot before the first optimisation step.
        * iteration N → snapshot after iteration N (1‑based).
        Files are written with the output prefix and zero‑padded four‑digit
        suffixes, e.g. for prefix "output_":
          output_tiles_0000.fits
          output_coverage_0000.png
          output_layers_0000.png
        """
        # Ensure the iteration number is an integer.
        it = int(iteration)
        # Zero‑pad to four digits.
        suffix = f"{it:04d}"
        prefix = self.output_prefix

        # ---- FITS checkpoint ----
        fits_name = f"{prefix}tiles_{suffix}.fits"
        # Use the existing save_checkpoint logic (iteration number stored in header).
        self.save_checkpoint(fits_name, it)

        # ---- Diagnostic PNGs (only when plotting is requested) ----
        # Write each plot directly to its suffixed filename (no intermediate rename).
        if self.plot:
            self.plot_diagnostics(out_path=f"{prefix}coverage_{suffix}.png")
            self._plot_allsky_layer0(out_path=f"{prefix}layers_{suffix}.png")


        
# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Optimize all-sky survey tiling coverage.")
    parser.add_argument("-c", "--camera", required=True, help="Camera FITS file")
    parser.add_argument("-t", "--tiles", required=True, help="Initial tile centres FITS file")
    parser.add_argument("-r", "--randoms", help="Random catalogue FITS file (optional)")
    parser.add_argument("--num-randoms", type=int, default=10_000_000, help="Number of randoms to generate")
    parser.add_argument("--max-drift", type=float, default=0.3, help="Max drift radius (deg)")
    parser.add_argument("--delta", type=float, default=0.1, help="Initial step size sigma (deg)")
    parser.add_argument("--shrink", type=float, help="Decay factor per iteration")
    parser.add_argument("--temp", type=float, default=0.1, help="Metropolis temperature")
    parser.add_argument("--iters", type=int, default=1000, help="Number of iterations")
    parser.add_argument("--workers", type=int, help="Number of parallel workers")
    parser.add_argument("--moves-per-tile", type=int, default=1,
                        help="Candidate moves proposed per tile per iteration; the best descent is kept (default 1)")
    parser.add_argument("--write-coverage", action="store_true", help="Write a 2nd HDU with coverage counts for each random")
    parser.add_argument("--plot", action="store_true", help="Generate diagnostic plots")
    parser.add_argument("--ra-center", type=float, help="RA center for zoom plot (deg)")
    parser.add_argument("--dec-center", type=float, help="Dec center for zoom plot (deg)")
    parser.add_argument("--diameter", type=float, help="Diameter for zoom plot (deg)")
    # New bounding‑box options to fix tiles outside a given range (override FIXED column)
    parser.add_argument("--fixed-ra-min", type=float, help="If set, tiles with RA < this value (deg) are fixed (outside the bound).")
    parser.add_argument("--fixed-ra-max", type=float, help="If set, tiles with RA > this value (deg) are fixed (outside the bound).")
    parser.add_argument("--fixed-dec-min", type=float, help="If set, tiles with Dec < this value (deg) are fixed (outside the bound).")
    parser.add_argument("--fixed-dec-max", type=float, help="If set, tiles with Dec > this value (deg) are fixed (outside the bound).")
    parser.add_argument("--force-fork", action="store_true", help="Force multiprocessing to use 'fork' (required for some macOS versions to support shared memory)")
    parser.add_argument("-o", "--output", default="output_",
                        help="Prefix for all output files (default 'output_'); e.g. "
                             "output_tiles_0000.fits, output_coverage_0000.png, "
                             "output_layers_0000.png, output_camera.png")

    args = parser.parse_args()

    # macOS fix for shared memory inheritance in multiprocessing
    if args.force_fork or sys.platform == 'darwin':
        try:
            mp.set_start_method('fork', force=True)
            print(f"Using 'fork' start method (detected platform: {sys.platform})")
        except RuntimeError as e:
            print(f"Warning: Could not set start method to fork: {e}")

    opt = TilingOptimizer(
        camera_fits=args.camera,
        tiles_fits=args.tiles,
        randoms_fits=args.randoms,
        num_randoms=args.num_randoms,
        max_drift=args.max_drift,
        delta=args.delta,
        shrink=args.shrink,
        temp=args.temp,
        iters=args.iters,
        workers=args.workers,
        moves_per_tile=args.moves_per_tile,
        output_prefix=args.output,
        write_random_coverage_hdu=args.write_coverage,
        plot=args.plot,
        ra_center=args.ra_center,
        dec_center=args.dec_center,
        diameter=args.diameter,
        fixed_ra_min=args.fixed_ra_min,
        fixed_ra_max=args.fixed_ra_max,
        fixed_dec_min=args.fixed_dec_min,
        fixed_dec_max=args.fixed_dec_max
    )

    opt.optimize()

if __name__ == "__main__":
    main()
