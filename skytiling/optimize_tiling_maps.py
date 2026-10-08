#!/usr/bin/env python3
"""
Tiling optimizer for cameras described by per-filter throughput maps.

A generalisation of ``optimize_tiling`` for cases like HSC with the MBQ1
quadrant filter: the camera is given by field-angle maps of relative
throughput for each sub-filter (``flats_to_throughput`` output), every tile
(pointing centre) is observed at a fixed set of instrument rotations, and the
coverage of a sky point in sub-filter f is

    c_f = sum over tiles, sum over rotations of  THRU_f(camera coords of the point)

i.e. the effective number of full-throughput exposures.  Each tile also
carries a rotation offset theta: it is observed at theta + each angle in
``--rotations`` (the set stays rigid, e.g. 90 deg apart, but can turn as a
whole; ``--delta-rot 0`` freezes it).  Tile centres and rotations are moved
by simulated annealing to minimise  sum_f Var(c_f)  over a catalogue of
random points.  The randoms can be restricted to a survey footprint, defined
as the union of discs around a list of pointing centres (e.g. the HSC-Niji
wide pointings); tiles are placed over the footprint plus a margin so its
edges are covered.

Typical use::

    optimize_tiling_maps -t data/subaru/hsc_mbq1_throughput.fits \\
        --footprint data/subaru/hsc_niji_wide_pointings.csv --footprint-radius 0.75 \\
        --ntiles 3208 --num-randoms 2000000 --iters 50 --plot -o out_hsc_

Initial tile positions are either a Fibonacci lattice clipped to the footprint
(+ margin), the pointing list itself with dither offsets (``--init-offsets``,
to start from or evaluate a dither pattern), or a tiles FITS file (``--tiles``,
e.g. a checkpoint to resume).  Each iteration writes
``<prefix>tiles_NNNN.fits`` (tileid, ra, dec, rot, ra0, dec0) and, with
``--plot``, coverage maps and histograms.
"""
import argparse
import os
import re
import sys
import time
import multiprocessing as mp

# One BLAS/OpenMP thread per worker process (see optimize_tiling.py).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from skytiling.fibonacci_tile import fibonacci_points, SPHERE_AREA_DEG2

DEG = np.pi / 180.0


# ----------------------------------------------------------------------------- geometry
def radec_to_xyz(ra, dec):
    ra = np.radians(ra); dec = np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def xyz_to_radec(xyz):
    ra = np.degrees(np.arctan2(xyz[..., 1], xyz[..., 0])) % 360.0
    dec = np.degrees(np.arcsin(np.clip(xyz[..., 2], -1, 1)))
    return ra, dec


def east_north(t):
    """Unit east and north vectors at unit vector t (array (...,3))."""
    t = np.asarray(t, float)
    z = np.array([0.0, 0.0, 1.0])
    east = np.cross(z, t); n = np.linalg.norm(east, axis=-1, keepdims=True)
    east = np.where(n > 1e-12, east / np.maximum(n, 1e-12), np.array([1.0, 0.0, 0.0]))
    north = np.cross(t, east)
    return east, north


def tangent_coords(pts, t, east, north):
    """Gnomonic (x=east, y=north) coordinates in degrees of unit vectors pts about centre t."""
    d = pts @ t
    d = np.where(d > 1e-6, d, 1e-6)
    return np.degrees((pts @ east) / d), np.degrees((pts @ north) / d)


# ----------------------------------------------------------------------------- camera
class ThroughputCamera:
    """Per-filter field-angle throughput maps with nearest-cell lookup."""

    def __init__(self, path, use_good=False):
        from astropy.io import fits
        with fits.open(path) as h:
            ph = h[0].header
            self.filters = [ph[f'FILT{i}'] for i in range(ph['NFILT'])]
            self.res = float(ph['RES']) / 60.0            # deg per cell
            self.half = float(ph['EXTENT']) / 60.0        # deg
            key = 'GOOD' if use_good else 'THRU'
            self.maps = np.stack([np.nan_to_num(h[f'{key}_{f}'].data.astype(np.float32), nan=0.0)
                                  for f in self.filters])   # (nf, N, N)
        self.nf = len(self.filters)
        self.n = self.maps.shape[1]
        any_map = self.maps.max(axis=0) > 0
        yy, xx = np.mgrid[0:self.n, 0:self.n]
        r = np.hypot((xx + 0.5) * self.res - self.half, (yy + 0.5) * self.res - self.half)
        self.max_radius = float(r[any_map].max())        # deg, outermost illuminated cell
        self.area = self.maps.sum(axis=(1, 2)) * self.res ** 2   # deg^2 per filter (throughput-weighted)

    def coverage(self, lx, ly, rotations_rad, theta=0.0):
        """Sum over rotations (each offset by theta) of the maps at camera coords;
        lx, ly in deg (tangent plane, x = east, y = north).  Returns (M, nf)."""
        out = np.zeros((lx.shape[0], self.nf), dtype=np.float32)
        for phi in rotations_rad:
            c, s = np.cos(phi + theta), np.sin(phi + theta)
            cx = lx * c + ly * s
            cy = -lx * s + ly * c
            ix = np.floor((cx + self.half) / self.res).astype(np.int64)
            iy = np.floor((cy + self.half) / self.res).astype(np.int64)
            ok = (ix >= 0) & (ix < self.n) & (iy >= 0) & (iy < self.n)
            out[ok] += self.maps[:, iy[ok], ix[ok]].T
        return out


# ----------------------------------------------------------------------------- footprint
class Footprint:
    """Union of discs of radius `radius_deg` around pointing centres."""

    def __init__(self, ra, dec, radius_deg):
        from scipy.spatial import cKDTree
        self.centers = radec_to_xyz(np.asarray(ra, float), np.asarray(dec, float))
        self.radius = radius_deg
        self.tree = cKDTree(self.centers)
        self.ra_range = (float(np.min(ra)), float(np.max(ra)))
        self.dec_range = (float(np.min(dec)), float(np.max(dec)))

    def contains(self, xyz, extra=0.0):
        chord = 2 * np.sin(np.radians(self.radius + extra) / 2)
        d, _ = self.tree.query(xyz, k=1)
        return d <= chord

    def interior(self, xyz, radius_deg=1.6, frac=0.9):
        """Points well inside the footprint: those with at least `frac` x the median
        number of pointing centres within `radius_deg` (edge points see fewer)."""
        chord = 2 * np.sin(np.radians(radius_deg) / 2)
        cnt = self.tree.query_ball_point(xyz, chord, return_length=True)
        return cnt >= frac * np.median(cnt)

    def random_points(self, n, rng, extra=0.0):
        """Uniform random points inside the footprint (+extra deg)."""
        dec_lo = max(-90.0, self.dec_range[0] - self.radius - extra - 0.1)
        dec_hi = min(90.0, self.dec_range[1] + self.radius + extra + 0.1)
        out = []
        got = 0
        while got < n:
            m = max(100000, int(1.5 * (n - got) / max(self._fill, 0.01)))
            ra = rng.uniform(0, 360, m)
            z = rng.uniform(np.sin(np.radians(dec_lo)), np.sin(np.radians(dec_hi)), m)
            dec = np.degrees(np.arcsin(z))
            xyz = radec_to_xyz(ra, dec)
            keep = self.contains(xyz, extra)
            self._fill = max(keep.mean(), 1e-4)
            out.append(xyz[keep]); got += keep.sum()
        return np.concatenate(out)[:n]

    _fill = 0.3

    def area(self, rng, nsamp=400000):
        """Monte-Carlo area in deg^2."""
        dec_lo = max(-90.0, self.dec_range[0] - self.radius - 0.1)
        dec_hi = min(90.0, self.dec_range[1] + self.radius + 0.1)
        ra = rng.uniform(0, 360, nsamp)
        z = rng.uniform(np.sin(np.radians(dec_lo)), np.sin(np.radians(dec_hi)), nsamp)
        band = 360.0 * np.degrees(np.sin(np.radians(dec_hi)) - np.sin(np.radians(dec_lo))) * (180 / np.pi) / 1.0
        band = 360.0 * (np.sin(np.radians(dec_hi)) - np.sin(np.radians(dec_lo))) * (180.0 / np.pi)
        frac = self.contains(radec_to_xyz(ra, np.degrees(np.arcsin(z)))).mean()
        return band * frac


class RectFootprint:
    """Union of RA/Dec rectangles [(ra_min, ra_max, dec_min, dec_max), ...] in deg
    (RA ranges may wrap through 0, e.g. 313.82..47.15).  Same interface as Footprint;
    `centers` is empty (use --init-lattice or --ntiles to seed tiles)."""

    def __init__(self, rects, names=None):
        self.rects = [tuple(float(v) for v in r) for r in rects]
        self.names = list(names) if names is not None else [f'rect{i}' for i in range(len(self.rects))]
        self.centers = np.zeros((0, 3))
        self.radius = 0.0
        self.dec_range = (min(r[2] for r in self.rects), max(r[3] for r in self.rects))
        self.ra_range = (0.0, 360.0)

    def _in_rect(self, ra, dec, rect, extra=0.0):
        ra_min, ra_max, dec_min, dec_max = rect
        # widen the RA bounds by extra / cos(dec) so the margin is a true angular distance
        cosd = np.maximum(np.cos(np.radians(np.clip(dec, -89.0, 89.0))), 1e-3)
        dra = (ra - ra_min) % 360.0
        width = (ra_max - ra_min) % 360.0
        in_ra = (dra <= width + extra / cosd) | (dra >= 360.0 - extra / cosd)
        return in_ra & (dec >= dec_min - extra) & (dec <= dec_max + extra)

    def contains(self, xyz, extra=0.0):
        ra, dec = xyz_to_radec(np.atleast_2d(xyz))
        out = np.zeros(len(ra), dtype=bool)
        for r in self.rects:
            out |= self._in_rect(ra, dec, r, extra)
        return out

    def interior(self, xyz, margin=0.75):
        """Points at least `margin` deg from the boundary of the union: the point and
        its 8 neighbours displaced by `margin` (N, S, E, W and diagonals) are all inside."""
        ra, dec = xyz_to_radec(np.atleast_2d(xyz))
        cosd = np.maximum(np.cos(np.radians(dec)), 1e-3)
        ok = self.contains(xyz)
        for ddec in (-margin, 0.0, margin):
            for dra in (-margin, 0.0, margin):
                if ddec == 0.0 and dra == 0.0:
                    continue
                ok &= self.contains(radec_to_xyz((ra + dra / cosd) % 360.0, np.clip(dec + ddec, -90, 90)))
        return ok

    _fill = 0.3

    def random_points(self, n, rng, extra=0.0):
        """Uniform random points inside the union of rectangles (+extra deg)."""
        dec_lo = max(-90.0, self.dec_range[0] - extra - 0.1)
        dec_hi = min(90.0, self.dec_range[1] + extra + 0.1)
        out = []
        got = 0
        while got < n:
            m = max(100000, int(1.5 * (n - got) / max(self._fill, 0.01)))
            ra = rng.uniform(0, 360, m)
            z = rng.uniform(np.sin(np.radians(dec_lo)), np.sin(np.radians(dec_hi)), m)
            dec = np.degrees(np.arcsin(z))
            xyz = radec_to_xyz(ra, dec)
            keep = self.contains(xyz, extra)
            self._fill = max(keep.mean(), 1e-4)
            out.append(xyz[keep]); got += keep.sum()
        return np.concatenate(out)[:n]

    def area(self, rng=None, nsamp=400000, extra=0.0):
        """Area of the union in deg^2 (Monte-Carlo, so overlapping rectangles are not double counted)."""
        rng = rng or np.random.default_rng(0)
        dec_lo = max(-90.0, self.dec_range[0] - extra - 0.1)
        dec_hi = min(90.0, self.dec_range[1] + extra + 0.1)
        ra = rng.uniform(0, 360, nsamp)
        z = rng.uniform(np.sin(np.radians(dec_lo)), np.sin(np.radians(dec_hi)), nsamp)
        band = 360.0 * (np.sin(np.radians(dec_hi)) - np.sin(np.radians(dec_lo))) * (180.0 / np.pi)
        return band * self.contains(radec_to_xyz(ra, np.degrees(np.arcsin(z))), extra).mean()

    def _bands(self):
        """Group rectangles whose Dec ranges touch or overlap into bands: [(dec_min, dec_max, [rects]), ...]."""
        def ra_overlap(r1, r2):
            a0, w1 = r1[0], (r1[1] - r1[0]) % 360.0
            b0, w2 = r2[0], (r2[1] - r2[0]) % 360.0
            return ((b0 - a0) % 360.0) <= w1 or ((a0 - b0) % 360.0) <= w2
        bands = []
        for r in sorted(self.rects, key=lambda r: r[2]):
            for b in bands:
                if r[2] <= b[1] + 1e-6 and r[3] >= b[0] - 1e-6 and any(ra_overlap(r, q) for q in b[2]):
                    b[0], b[1] = min(b[0], r[2]), max(b[1], r[3]); b[2].append(r); break
            else:
                bands.append([r[2], r[3], [r]])
        return bands

    def lattice(self, dra, ddec, extra=0.0):
        """Pointing centres on rows of constant Dec.  Rectangles that touch in Dec form one
        band; each band gets round(height / ddec) rows, equally spaced and centred in the
        band (the first row half a spacing above dec_min).  Along a row, the RA extent is
        the union of the rectangles containing that Dec and the step is dra / cos(dec)
        adjusted so a whole number of points is centred in the extent; alternate rows are
        offset by half a step.  Returns xyz (n, 3), those inside the union (+extra)."""
        pts = []
        for dec_lo, dec_hi, rects in self._bands():
            nrow = max(1, int(round((dec_hi - dec_lo) / ddec)))
            spacing = (dec_hi - dec_lo) / nrow
            for k in range(nrow):
                dec = dec_lo + (k + 0.5) * spacing
                # RA intervals of the rectangles containing this row, merged (unwrapped from the first ra_min)
                ivs = []
                for r in rects:
                    if r[2] - 1e-6 <= dec <= r[3] + 1e-6:
                        a = (r[0] - rects[0][0] + 180.0) % 360.0 - 180.0     # start, unwrapped about the first rect
                        ivs.append((a, a + (r[1] - r[0]) % 360.0))
                ivs.sort(); merged = []
                for a, b in ivs:
                    if merged and a <= merged[-1][1] + 1e-6:
                        merged[-1][1] = max(merged[-1][1], b)
                    else:
                        merged.append([a, b])
                step0 = dra / max(np.cos(np.radians(dec)), 1e-3)
                for a, b in merged:
                    ncol = max(1, int(round((b - a) / step0)))
                    step = (b - a) / ncol
                    for j in range(ncol):
                        off = (j + 0.5 + (0.5 if k % 2 else 0.0)) * step
                        if off < b - a:
                            pts.append(((rects[0][0] + a + off) % 360.0, dec))
        pts = np.unique(np.round(np.array(pts), 6), axis=0)
        xyz = radec_to_xyz(pts[:, 0], pts[:, 1])
        return xyz[self.contains(xyz, extra)]


def parse_rects(spec):
    """'ra_min,ra_max,dec_min,dec_max[;...]' or a CSV file with those columns (and an optional name)."""
    if os.path.exists(spec):
        import csv
        rows = list(csv.DictReader(l for l in open(spec) if not l.startswith('#')))
        rects = [(r['ra_min'], r['ra_max'], r['dec_min'], r['dec_max']) for r in rows]
        names = [r.get('name', f'rect{i}') for i, r in enumerate(rows)]
    elif spec.endswith('.csv') or '/' in spec:
        sys.exit(f'--rects: file not found: {spec} (relative to {os.getcwd()})')
    else:
        rects = [tuple(v for v in part.split(',')) for part in spec.split(';') if part.strip()]
        names = None
    return RectFootprint(rects, names)


# ----------------------------------------------------------------------------- objectives
def energy_from_sums(sum_c, sum_c2, n, objective='var', power=1.0):
    """Energy of a coverage distribution from its per-filter sums of c and c^2 over n randoms.
    'var': sum_f (S2 - S1^2/n) = n sum_f Var(c_f) (the original objective; scale ~ n).
    'cv' : sum_f rms_f / mean_f^power (scale ~ 1; rewards a higher mean inside the footprint,
           i.e. penalises coverage spilled outside it)."""
    sum_c = np.asarray(sum_c, dtype=np.float64); sum_c2 = np.asarray(sum_c2, dtype=np.float64)
    if objective == 'var':
        return float(np.sum(sum_c2 - sum_c ** 2 / n))
    mean = np.maximum(sum_c / n, 1e-12)
    var = np.maximum(sum_c2 / n - mean ** 2, 0.0)
    return float(np.sum(np.sqrt(var) / mean ** power))


# ----------------------------------------------------------------------------- optimizer
# Module-level state inherited by forked worker processes (never pickled per task):
# shared randoms, shared coverage, candidate index lists, camera, initial tile positions.
_G = {}


def _tile_coverage(cam, rot, t, pts, theta=0.0):
    east, north = east_north(t)
    lx, ly = tangent_coords(pts, t, east, north)
    return cam.coverage(lx, ly, rot, theta)


def _worker(task):
    """Process a chunk of tiles from one layer: propose moves, accept, update shared coverage.

    task = (items, sum_c, sum_c2, sigma, sigma_rot, seed) with items = [(i, tile_xyz, theta), ...].
    Tiles in one layer have disjoint candidate sets, so workers never write the
    same coverage rows.  Returns [(i, new_xyz, new_theta, accepted, d_sum, d_sum2, move_deg, rot_deg), ...].
    """
    items, sum_c, sum_c2, sigma, sigma_rot, seed = task
    cam, rot = _G['cam'], _G['rot']
    n_rand, nf = _G['n_rand'], cam.nf
    cov = np.frombuffer(_G['cov_raw'], dtype=np.float32).reshape(n_rand, nf)
    rand = np.frombuffer(_G['xyz_raw'], dtype=np.float64).reshape(n_rand, 3)
    cand, tiles0 = _G['cand'], _G['tiles0']
    max_drift, temp, k = _G['max_drift'], _G['temp'], _G['moves_per_tile']
    objective, power = _G.get('objective', 'var'), _G.get('cv_power', 1.0)
    keep_fp, keep_margin = _G.get('keep_fp'), _G.get('keep_margin', 0.0)
    rng = np.random.default_rng(seed)
    sum_c = np.array(sum_c, dtype=np.float64); sum_c2 = np.array(sum_c2, dtype=np.float64)
    e_cur = energy_from_sums(sum_c, sum_c2, n_rand, objective, power) if objective != 'var' else 0.0
    out = []
    for i, t, th in items:
        c = cand[i]
        if len(c) == 0:
            out.append((i, t, th, False, None, None, 0.0, 0.0)); continue
        pts = rand[c]
        cur = cov[c].astype(np.float64)
        old = _tile_coverage(cam, rot, t, pts, th).astype(np.float64)
        east, north = east_north(t)
        best = None
        for _ in range(k):
            dx, dy = rng.normal(0, sigma, 2) if sigma > 0 else (0.0, 0.0)
            p = t + dx * east + dy * north
            p /= np.linalg.norm(p)
            d0 = np.arccos(np.clip(p @ tiles0[i], -1, 1))
            if d0 > max_drift:
                p = tiles0[i] + (p - tiles0[i]) * (max_drift / d0)
                p /= np.linalg.norm(p)
            if keep_fp is not None and not keep_fp.contains(p[None, :], keep_margin)[0]:
                continue                                   # proposal would leave footprint + margin
            th_new = th + (rng.normal(0, sigma_rot) if sigma_rot > 0 else 0.0)
            new = _tile_coverage(cam, rot, p, pts, th_new).astype(np.float64)
            delta = new - old
            d_sum = delta.sum(axis=0)
            d_sum2 = (2 * cur * delta + delta ** 2).sum(axis=0)
            if objective == 'var':
                d_e = float(np.sum(d_sum2 - (2 * sum_c * d_sum + d_sum ** 2) / n_rand))
            else:
                d_e = energy_from_sums(sum_c + d_sum, sum_c2 + d_sum2, n_rand, objective, power) - e_cur
            if best is None or d_e < best[0]:
                best = (d_e, p, th_new, delta, d_sum, d_sum2)
        if best is None:
            out.append((i, t, th, False, None, None, 0.0, 0.0)); continue
        d_e, p, th_new, delta, d_sum, d_sum2 = best
        tscale = 1e6 if objective == 'var' else 1.0
        accept = d_e < 0 or (temp > 0 and rng.random() < np.exp(-d_e / (temp * tscale)))
        if accept:
            cov[c] += delta.astype(np.float32)
            sum_c += d_sum; sum_c2 += d_sum2
            if objective != 'var':
                e_cur = energy_from_sums(sum_c, sum_c2, n_rand, objective, power)
            out.append((i, p, th_new, True, d_sum, d_sum2,
                        float(np.degrees(np.arccos(np.clip(p @ t, -1, 1)))), float(np.degrees(th_new - th))))
        else:
            out.append((i, t, th, False, None, None, 0.0, 0.0))
    return out


class MapTilingOptimizer:
    def __init__(self, camera, randoms_xyz, tiles_xyz, rotations_deg, max_drift, delta,
                 temp, moves_per_tile, seed=0, workers=1, theta_deg=None, delta_rot=0.0,
                 objective='var', cv_power=1.0, keep_fp=None, keep_margin=0.0):
        from scipy.spatial import cKDTree
        self.cam = camera
        self.objective, self.cv_power = objective, cv_power
        self.keep_fp, self.keep_margin = keep_fp, keep_margin
        self.rot = [np.radians(r) for r in rotations_deg]
        self.n_rand = len(randoms_xyz)
        self.tiles = np.ascontiguousarray(tiles_xyz, dtype=np.float64)
        self.tiles0 = self.tiles.copy()
        self.n_tiles = len(self.tiles)
        self.theta = np.zeros(self.n_tiles) if theta_deg is None else np.radians(np.asarray(theta_deg, float))
        self.max_drift = np.radians(max_drift)
        self.sigma = np.radians(delta)
        self.sigma_rot = np.radians(delta_rot)
        self.temp = temp
        self.moves_per_tile = max(1, int(moves_per_tile))
        self.workers = max(1, int(workers))
        self.rng = np.random.default_rng(seed)
        self.interior = np.ones(self.n_rand, dtype=bool)
        self.pool = None

        # shared randoms and coverage (RawArrays inherited by forked workers)
        self._xyz_raw = mp.RawArray('d', self.n_rand * 3)
        self.rand = np.frombuffer(self._xyz_raw, dtype=np.float64).reshape(self.n_rand, 3)
        self.rand[:] = randoms_xyz
        self._cov_raw = mp.RawArray('f', self.n_rand * self.cam.nf)
        self.cov = np.frombuffer(self._cov_raw, dtype=np.float32).reshape(self.n_rand, self.cam.nf)
        self.cov[:] = 0

        t0 = time.time()
        self.tree = cKDTree(self.rand)
        chord = 2 * np.sin((np.radians(self.cam.max_radius) + self.max_drift) / 2)
        self.cand = [np.array(c, dtype=np.int64) for c in self.tree.query_ball_point(self.tiles, chord)]
        print(f'[timing] candidate sets in {time.time() - t0:.1f} s; '
              f'mean {np.mean([len(c) for c in self.cand]):.0f} randoms per tile', flush=True)

        t0 = time.time()
        for i in range(self.n_tiles):
            c = self.cand[i]
            if len(c):
                self.cov[c] += self._tile_coverage(self.tiles[i], self.rand[c], self.theta[i])
        print(f'[timing] initial coverage in {time.time() - t0:.1f} s', flush=True)
        self._refresh_sums()

        t0 = time.time()
        self.layers = self._build_layers()
        print(f'[timing] {len(self.layers)} layers (mean {self.n_tiles / len(self.layers):.1f} tiles) in {time.time() - t0:.1f} s', flush=True)

    def _tile_coverage(self, t, pts, theta=0.0):
        return _tile_coverage(self.cam, self.rot, t, pts, theta)

    def _build_layers(self):
        """Greedy colouring so that no two tiles closer than 2*(max_radius + max_drift)
        share a layer: their candidate sets are then disjoint and can be updated in parallel."""
        from scipy.spatial import cKDTree
        sep = 2 * (np.radians(self.cam.max_radius) + self.max_drift)
        chord = 2.0 if sep >= np.pi else 2 * np.sin(sep / 2)
        tree = cKDTree(self.tiles0)
        nb = tree.query_ball_point(self.tiles0, chord)
        order = sorted(range(self.n_tiles), key=lambda i: len(self.cand[i]), reverse=True)
        colour = np.full(self.n_tiles, -1)
        for i in order:
            used = {colour[j] for j in nb[i] if colour[j] >= 0}
            c = 0
            while c in used:
                c += 1
            colour[i] = c
        layers = [[] for _ in range(colour.max() + 1)]
        for i in order:
            layers[colour[i]].append(i)
        return layers

    def _start_pool(self):
        """Fork the worker pool with the shared state in place (call after the
        randoms, coverage, candidates and tiles0 are set and will not be replaced)."""
        if self.workers <= 1:
            return
        _G.update(cam=self.cam, rot=self.rot, n_rand=self.n_rand, cov_raw=self._cov_raw, xyz_raw=self._xyz_raw,
                  cand=self.cand, tiles0=self.tiles0, max_drift=self.max_drift, temp=self.temp,
                  moves_per_tile=self.moves_per_tile, objective=self.objective, cv_power=self.cv_power,
                  keep_fp=self.keep_fp, keep_margin=self.keep_margin)
        ctx = mp.get_context('fork')
        self.pool = ctx.Pool(self.workers)

    def _refresh_sums(self):
        self.sum_c = self.cov.sum(axis=0, dtype=np.float64)
        self.sum_c2 = (self.cov.astype(np.float64) ** 2).sum(axis=0)

    def stats(self):
        mean = self.sum_c / self.n_rand
        var = np.maximum(self.sum_c2 / self.n_rand - mean ** 2, 0)
        return mean, np.sqrt(var)

    def energy(self):
        return energy_from_sums(self.sum_c, self.sum_c2, self.n_rand, self.objective, self.cv_power)

    def iterate(self):
        """One sweep over all layers; returns summary dict."""
        if self.workers > 1 and self.pool is None:
            self._start_pool()
        _G.update(cam=self.cam, rot=self.rot, n_rand=self.n_rand, cov_raw=self._cov_raw, xyz_raw=self._xyz_raw,
                  cand=self.cand, tiles0=self.tiles0, max_drift=self.max_drift, temp=self.temp,
                  moves_per_tile=self.moves_per_tile, objective=self.objective, cv_power=self.cv_power,
                  keep_fp=self.keep_fp, keep_margin=self.keep_margin)
        acc = imp = rej = 0
        dists, rots = [], []
        for layer in self.layers:
            layer = [int(i) for i in self.rng.permutation(layer)]
            nchunk = min(self.workers, len(layer)) if self.workers > 1 else 1
            chunks = [layer[j::nchunk] for j in range(nchunk)]
            tasks = [([(i, self.tiles[i], float(self.theta[i])) for i in ch], self.sum_c.copy(), self.sum_c2.copy(),
                      self.sigma, self.sigma_rot, int(self.rng.integers(2 ** 31))) for ch in chunks if ch]
            results = self.pool.map(_worker, tasks) if self.pool else [_worker(t) for t in tasks]
            for res in results:
                for i, p, th, accepted, d_sum, d_sum2, move, drot in res:
                    if accepted:
                        self.tiles[i] = p; self.theta[i] = th
                        self.sum_c += d_sum; self.sum_c2 += d_sum2
                        dists.append(move); rots.append(abs(drot)); acc += 1
                    else:
                        rej += 1
        self._refresh_sums()   # exact sums (guards float32 accumulation and stale sums in workers)
        mean, rms = self.stats()
        return dict(accepted=acc, improved=acc, rejected=rej, mean=mean, rms=rms,
                    move=(np.min(dists), np.median(dists), np.max(dists)) if dists else (0, 0, 0),
                    rot=(np.median(rots), np.max(rots)) if rots else (0, 0))

    def close(self):
        if self.pool:
            self.pool.close(); self.pool.join(); self.pool = None

    # ----------------------------------------------------------------- output
    def write_tiles(self, path, iteration):
        from astropy.io import fits
        ra, dec = xyz_to_radec(self.tiles); ra0, dec0 = xyz_to_radec(self.tiles0)
        cols = [fits.Column(name='tileid', format='K', array=np.arange(self.n_tiles) + 1),
                fits.Column(name='ra', format='D', array=ra), fits.Column(name='dec', format='D', array=dec),
                fits.Column(name='rot', format='D', unit='deg', array=np.degrees(self.theta) % 360.0),
                fits.Column(name='ra0', format='D', array=ra0), fits.Column(name='dec0', format='D', array=dec0)]
        if getattr(self, 'pid', None) is not None:
            cols.append(fits.Column(name='pid', format='K', array=np.asarray(self.pid, dtype=np.int64)))
        if getattr(self, 'pname', None) is not None:
            w = max(len(str(v)) for v in self.pname)
            cols.append(fits.Column(name='pname', format=f'{w}A', array=np.asarray(self.pname, dtype=f'U{w}')))
        hdu = fits.BinTableHDU.from_columns(cols, name='TILES')
        hdu.header['ITER'] = iteration
        hdu.header['OBJECTIV'] = (self.objective + (f' p={self.cv_power:g}' if self.objective == 'cv' else ''), 'annealing objective')
        hdu.header['ENERGY'] = (self.energy(), 'objective value')
        hdu.header['ROTS'] = (','.join(f'{np.degrees(r):g}' for r in self.rot), 'rotation set added to rot [deg]')
        mean, rms = self.stats()
        for k, f in enumerate(self.cam.filters):
            hdu.header[f'MEAN_{f}'] = (float(mean[k]), 'mean coverage'); hdu.header[f'RMS_{f}'] = (float(rms[k]), 'rms coverage')
        hdu.writeto(path, overwrite=True)

    def report(self, label=''):
        for name, sel in (('all randoms', slice(None)), (f'interior ({100 * self.interior.mean():.0f}%)', self.interior)):
            c = self.cov[sel]
            mean = c.mean(axis=0); rms = c.std(axis=0)
            print(f'{label} {name}: ' + '; '.join(
                f'{f} {mean[k]:.3f}+-{rms[k]:.3f} ({rms[k] / mean[k]:.3f}) <1.5 {100 * np.mean(c[:, k] < 1.5):.2f}% <0.5 {100 * np.mean(c[:, k] < 0.5):.3f}%'
                for k, f in enumerate(self.cam.filters)))

    def grid_coverage(self, ra_center, dec_center, diameter, res_arcmin=0.5):
        """Coverage (n_pix, n_pix, nf) evaluated on a tangent-plane grid about (ra_center, dec_center)."""
        n = int(round(diameter * 60 / res_arcmin))
        g = (np.arange(n) + 0.5) * res_arcmin / 60.0 - diameter / 2
        X, Y = np.meshgrid(g, g)                       # X = east offset [deg], Y = north
        c0 = radec_to_xyz(ra_center, dec_center)
        e0, n0 = east_north(c0)
        # inverse gnomonic projection of the grid
        xr, yr = np.radians(X.ravel()), np.radians(Y.ravel())
        pts = c0[None, :] + xr[:, None] * e0[None, :] + yr[:, None] * n0[None, :]
        pts /= np.linalg.norm(pts, axis=1, keepdims=True)
        cov = np.zeros((pts.shape[0], self.cam.nf), dtype=np.float32)
        near = np.arccos(np.clip(self.tiles @ c0, -1, 1)) < np.radians(diameter / 2 * 1.5 + self.cam.max_radius)
        for t, th in zip(self.tiles[near], self.theta[near]):
            sel = (pts @ t) > np.cos(np.radians(self.cam.max_radius))
            if sel.any():
                cov[sel] += self._tile_coverage(t, pts[sel], th)
        return X, Y, cov.reshape(n, n, self.cam.nf)

    def plot(self, path, ra_center, dec_center, diameter, res_arcmin=0.5):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        X, Y, cov = self.grid_coverage(ra_center, dec_center, diameter, res_arcmin)
        h = diameter / 2
        nf = self.cam.nf
        fig, axes = plt.subplots(2, nf, figsize=(4.6 * nf, 8.5))
        vmax = float(np.percentile(self.cov, 99.5)) * 1.05
        tra, tdec = xyz_to_radec(self.tiles)
        tdra = ((tra - ra_center + 180) % 360 - 180) * np.cos(np.radians(dec_center)); tddec = tdec - dec_center
        tsel = (np.abs(tdra) < h) & (np.abs(tddec) < h)
        mean, rms = self.stats()
        for k, f in enumerate(self.cam.filters):
            ax = axes[0, k]
            im = ax.imshow(cov[:, :, k], origin='lower', extent=[h, -h, -h, h], vmin=0, vmax=vmax, cmap='viridis')
            ax.plot(tdra[tsel], tddec[tsel], 'w+', ms=4, mew=0.8)
            ax.set_title(f'{f}: coverage around ({ra_center}, {dec_center})'); ax.set_xlabel('dRA cos(dec) [deg]'); ax.set_ylabel('dDec [deg]')
            plt.colorbar(im, ax=ax, fraction=0.046)
            ax = axes[1, k]
            ax.hist(self.cov[:, k], bins=np.arange(0, vmax + 0.25, 0.25), color=f'C{k}')
            ax.set_title(f'{f}: mean {mean[k]:.2f}, rms {rms[k]:.2f}, rms/mean {rms[k] / mean[k]:.3f}')
            ax.set_xlabel('coverage (sum of throughput)'); ax.set_ylabel('randoms')
        plt.tight_layout(); plt.savefig(path, dpi=80); plt.close(fig)


# ----------------------------------------------------------------------------- driver
def parse_resume_iteration(path):
    m = re.search(r'tiles_(\d+)\.fits$', os.path.basename(path))
    return int(m.group(1)) if m else None


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-t', '--throughput', required=True, help='flats_to_throughput FITS file')
    p.add_argument('--use-good', action='store_true', help='use the 0/1 GOOD maps instead of THRU')
    p.add_argument('--rotations', default='0,90,180,270', help='instrument rotations per tile [deg]')
    p.add_argument('--footprint', help='CSV with ra,dec columns of pointing centres defining the survey footprint')
    p.add_argument('--footprint-radius', type=float, default=0.75, help='disc radius around each centre [deg]')
    p.add_argument('--rects', help='rectangular footprint instead: "ra_min,ra_max,dec_min,dec_max;..." [deg] or a CSV '
                                   'with those columns (+ name); RA ranges may wrap through 0')
    p.add_argument('--interior-margin', type=float, default=0.75,
                   help='for --rects: "interior" statistics exclude points within this of the boundary [deg]')
    p.add_argument('--margin', type=float, default=0.4, help='tiles are placed out to this far beyond the footprint [deg]')
    p.add_argument('--randoms', help='FITS file of randoms (RA, DEC) to reuse')
    p.add_argument('--num-randoms', type=int, default=2_000_000)
    p.add_argument('--tiles', help='initial tiles FITS (ra, dec); a <prefix>tiles_NNNN.fits checkpoint resumes')
    p.add_argument('--ntiles', type=int, help='number of Fibonacci tiles over footprint+margin')
    p.add_argument('--tilearea', type=float, help='Fibonacci tile area [deg^2] (alternative to --ntiles)')
    p.add_argument('--init-offsets', help='start from the footprint pointings with these dither offsets '
                                           '"dx,dy dx,dy ..." [arcmin, +x = +RA]; e.g. Atsushi phase 1')
    p.add_argument('--init-lattice', help='for --rects: pointing centres on rows of constant Dec, "dra,ddec" [deg] '
                                          '(e.g. "1.299,1.125" = the HSC-SSP Wide lattice); combine with --init-offsets')
    p.add_argument('--init-centers', help='CSV of pointing centres (ra, dec; optional id) to start from, keeping those inside '
                                          'the footprint (+ --margin); pid = id column or row index; combine with --init-offsets')
    p.add_argument('--objective', choices=['var', 'cv'], default='var',
                   help='var: sum_f N Var(c_f) (default); cv: sum_f rms_f / mean_f^power, which also rewards '
                        'coverage kept inside the footprint')
    p.add_argument('--cv-power', type=float, default=1.0, help='power of the mean in the cv objective')
    p.add_argument('--keep-inside', type=float, help='reject moves that take a tile centre further than this '
                                                     'beyond the footprint [deg]')
    p.add_argument('--max-drift', type=float, default=1.0, help='max drift from initial position [deg]')
    p.add_argument('--delta', type=float, default=0.1, help='position proposal sigma [deg]')
    p.add_argument('--delta-rot', type=float, default=0.0,
                   help='rotation-offset proposal sigma [deg]; 0 keeps every tile at the --rotations set')
    p.add_argument('--init-rot', default='zero', help='initial rotation offsets: "zero", "random", or a value in deg')
    p.add_argument('--shrink', type=float, default=1.0, help='multiply delta by this each iteration')
    p.add_argument('--temp', type=float, default=0.0, help='Metropolis temperature')
    p.add_argument('--moves-per-tile', type=int, default=3)
    p.add_argument('--iters', type=int, default=20)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--workers', type=int, default=0,
                   help='worker processes (0 = NCPUS / SLURM_CPUS_PER_TASK / all CPUs)')
    p.add_argument('-o', '--output', default='output_', help='output prefix')
    p.add_argument('--plot', action='store_true')
    p.add_argument('--ra-center', type=float); p.add_argument('--dec-center', type=float)
    p.add_argument('--diameter', type=float, default=5.0, help='zoom plot size [deg]')
    return p.parse_args(argv)


def main(argv=None):
    from astropy.io import fits
    opts = parse_arguments(argv)
    rng = np.random.default_rng(opts.seed)
    cam = ThroughputCamera(opts.throughput, use_good=opts.use_good)
    rot = [float(v) for v in opts.rotations.split(',')]
    print(f'camera: filters {cam.filters}, max radius {cam.max_radius:.3f} deg, '
          f'area per filter {np.round(cam.area, 4)} deg^2 ({"GOOD" if opts.use_good else "THRU"})')

    # footprint and randoms
    fp = None
    if opts.rects:
        fp = parse_rects(opts.rects)
        area = fp.area(rng)
        print(f'footprint: {len(fp.rects)} rectangle(s) {fp.names}, area {area:.1f} deg^2')
    elif opts.footprint:
        import csv
        rows = list(csv.DictReader(open(opts.footprint)))
        fra = np.array([float(r['ra']) for r in rows]); fdec = np.array([float(r['dec']) for r in rows])
        fp = Footprint(fra, fdec, opts.footprint_radius)
        area = fp.area(rng)
        print(f'footprint: {len(rows)} centres, radius {opts.footprint_radius} deg, area {area:.1f} deg^2')
    if opts.randoms:
        d = fits.getdata(opts.randoms)
        rand = radec_to_xyz(d['RA'], d['DEC'])
        print(f'read {len(rand)} randoms from {opts.randoms}')
    else:
        if fp is None:
            sys.exit('need --footprint (or --randoms) to define where coverage is evaluated')
        rand = fp.random_points(opts.num_randoms, rng)
        ra, dec = xyz_to_radec(rand)
        fits.BinTableHDU.from_columns([fits.Column(name='RA', format='D', array=ra),
                                       fits.Column(name='DEC', format='D', array=dec)]).writeto(
            f'{opts.output}randoms_{len(rand)}.fits', overwrite=True)
        print(f'generated {len(rand)} randoms in the footprint -> {opts.output}randoms_{len(rand)}.fits')

    # initial tiles
    start_iter = 0
    theta = None
    pid = None          # pointing id (index of the undithered centre) when tiles come from centres x offsets
    pname = None        # pointing name from --init-centers, if its CSV has one
    if opts.tiles:
        d = fits.getdata(opts.tiles)
        tiles = radec_to_xyz(d['ra'], d['dec'])
        if 'rot' in d.names:
            theta = np.array(d['rot'], float)
        if 'pid' in d.names:
            pid = np.array(d['pid'], int)
        if 'pname' in d.names:
            pname = np.array(d['pname'])
        it = parse_resume_iteration(opts.tiles)
        if it is not None and 'ra0' in d.names:
            start_iter = it
            tiles0 = radec_to_xyz(d['ra0'], d['dec0'])
        else:
            tiles0 = tiles.copy()
        print(f'read {len(tiles)} tiles from {opts.tiles}' + (f', resuming at iteration {it}' if it is not None else ''))
    elif opts.init_offsets or opts.init_lattice or opts.init_centers:
        center_ids = None
        if opts.init_centers:
            import csv
            crows = list(csv.DictReader(l for l in open(opts.init_centers) if not l.startswith('#')))
            cxyz = radec_to_xyz(np.array([float(r['ra']) for r in crows]), np.array([float(r['dec']) for r in crows]))
            cid = np.array([int(r['id']) if 'id' in r else i for i, r in enumerate(crows)])
            keep = fp.contains(cxyz, extra=opts.margin) if fp is not None else np.ones(len(cxyz), bool)
            centers, center_ids = cxyz[keep], cid[keep]
            if 'name' in crows[0]:
                pname = np.repeat(np.array([r['name'] for r in crows])[keep], max(1, len(opts.init_offsets.split()) if opts.init_offsets else 1))
            print(f'{len(centers)} of {len(crows)} centres from {opts.init_centers} inside footprint + {opts.margin} deg')
        elif opts.init_lattice:
            if not isinstance(fp, RectFootprint):
                sys.exit('--init-lattice needs --rects')
            dra, ddec = (float(v) for v in opts.init_lattice.split(','))
            centers = fp.lattice(dra, ddec, extra=opts.margin)
            print(f'lattice of {len(centers)} pointing centres ({dra} x {ddec} deg) inside footprint + {opts.margin} deg')
        else:
            centers = fp.centers
        off = np.radians(np.array([[float(v) for v in o.split(',')] for o in opts.init_offsets.split()]) / 60.0) \
            if opts.init_offsets else np.zeros((1, 2))
        east, north = east_north(centers)
        tiles = (centers[:, None, :] + off[None, :, 0, None] * east[:, None, :] + off[None, :, 1, None] * north[:, None, :]).reshape(-1, 3)
        tiles /= np.linalg.norm(tiles, axis=1, keepdims=True)
        tiles0 = tiles.copy()
        pid = np.repeat(np.arange(len(centers)) if center_ids is None else center_ids, len(off))
        print(f'{len(tiles)} tiles = {len(centers)} pointings x {len(off)} offsets')
    else:
        if fp is None:
            sys.exit('need --footprint for a Fibonacci initial tiling')
        rot0 = rng.uniform(0, 360)   # lattice rotation about the pole, so its seam is not aligned with the footprint

        def clipped_lattice(tilearea):
            n_all = int(round(SPHERE_AREA_DEG2 / tilearea))
            _, fra_all, fdec_all = fibonacci_points(n_all)
            xyz = radec_to_xyz((fra_all + rot0) % 360, fdec_all)
            return n_all, xyz[fp.contains(xyz, extra=opts.margin)]

        if opts.ntiles:
            # bisect the lattice density until exactly --ntiles fall inside footprint + margin
            area_ext = fp.area(rng, extra=opts.margin) if isinstance(fp, RectFootprint) \
                else Footprint(fra, fdec, opts.footprint_radius + opts.margin).area(rng)
            lo, hi = 0.5 * area_ext / opts.ntiles, 2.0 * area_ext / opts.ntiles
            for _ in range(60):
                tilearea = 0.5 * (lo + hi)
                n_all, tiles = clipped_lattice(tilearea)
                if len(tiles) == opts.ntiles:
                    break
                if len(tiles) > opts.ntiles:
                    lo = tilearea
                else:
                    hi = tilearea
        elif opts.tilearea:
            tilearea = opts.tilearea
            n_all, tiles = clipped_lattice(tilearea)
        else:
            sys.exit('give --ntiles or --tilearea for the Fibonacci initial tiling')
        tiles0 = tiles.copy()
        print(f'Fibonacci lattice: {tilearea:.4f} deg^2 per tile ({n_all} all-sky), {len(tiles)} inside footprint + {opts.margin} deg margin')
    print(f'expected mean coverage per filter ~ {np.round(len(tiles) * len(rot) * cam.area / (fp.area(rng) if fp else 1), 3)} '
          f'(tiles x rotations x area / footprint area; edge tiles cover outside)')

    if theta is None:
        if opts.init_rot == 'random':
            theta = rng.uniform(0, 360, len(tiles))
        elif opts.init_rot == 'zero':
            theta = np.zeros(len(tiles))
        else:
            theta = np.full(len(tiles), float(opts.init_rot))
    workers = opts.workers or int(os.environ.get('NCPUS', os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1)))
    print(f'using {workers} worker process(es)')
    opt = MapTilingOptimizer(cam, rand, tiles, rot, opts.max_drift, opts.delta, opts.temp, opts.moves_per_tile,
                             seed=opts.seed, workers=workers, theta_deg=theta, delta_rot=opts.delta_rot,
                             objective=opts.objective, cv_power=opts.cv_power,
                             keep_fp=fp if opts.keep_inside is not None else None, keep_margin=opts.keep_inside or 0.0)
    opt.tiles0[:] = tiles0
    opt.pid = pid; opt.pname = pname
    if fp is not None:
        opt.interior = fp.interior(rand, opts.interior_margin) if isinstance(fp, RectFootprint) else fp.interior(rand)
    print(f'objective: {opts.objective}' + (f' (power {opts.cv_power})' if opts.objective == 'cv' else '')
          + (f'; tile centres kept within {opts.keep_inside} deg of the footprint' if opts.keep_inside is not None else ''))
    opt.report('initial')
    if opts.ra_center is None:
        rc, dc = xyz_to_radec(np.mean(opt.tiles, axis=0)[None])
        opts.ra_center, opts.dec_center = round(float(rc[0]), 2), round(float(dc[0]), 2)
    if start_iter == 0:
        opt.write_tiles(f'{opts.output}tiles_0000.fits', 0)
        if opts.plot:
            opt.plot(f'{opts.output}coverage_0000.png', opts.ra_center, opts.dec_center, opts.diameter)
    for it in range(start_iter + 1, start_iter + opts.iters + 1):
        t0 = time.time()
        r = opt.iterate()
        s = '; '.join(f'{f} {r["mean"][k]:.3f}+-{r["rms"][k]:.3f}' for k, f in enumerate(cam.filters))
        print(f'Iter {it:4d} | {s} | sum rms/mean {np.sum(r["rms"] / r["mean"]):.4f} | E {opt.energy():.5g} | sigma {np.degrees(opt.sigma):.4f} '
              f'| acc {r["accepted"]} rej {r["rejected"]} | move [{r["move"][0]:.4f}, {r["move"][1]:.4f}, {r["move"][2]:.4f}] '
              f'| rot [{r["rot"][0]:.2f}, {r["rot"][1]:.2f}] | {time.time() - t0:.1f} s', flush=True)
        opt.write_tiles(f'{opts.output}tiles_{it:04d}.fits', it)
        if opts.plot:
            opt.plot(f'{opts.output}coverage_{it:04d}.png', opts.ra_center, opts.dec_center, opts.diameter)
        opt.sigma *= opts.shrink; opt.sigma_rot *= opts.shrink
    opt.close()
    opt.report('final')


if __name__ == '__main__':
    main()
