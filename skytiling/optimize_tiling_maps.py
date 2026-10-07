#!/usr/bin/env python3
"""
Tiling optimizer for cameras described by per-filter throughput maps.

A generalisation of ``optimize_tiling`` for cases like HSC with the MBQ1
quadrant filter: the camera is given by field-angle maps of relative
throughput for each sub-filter (``flats_to_throughput`` output), every tile
(pointing centre) is observed at a fixed set of instrument rotations, and the
coverage of a sky point in sub-filter f is

    c_f = sum over tiles, sum over rotations of  THRU_f(camera coords of the point)

i.e. the effective number of full-throughput exposures.  Tile centres are
moved by simulated annealing to minimise  sum_f Var(c_f)  over a catalogue of
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
``<prefix>tiles_NNNN.fits`` (tileid, ra, dec, ra0, dec0) and, with
``--plot``, coverage maps and histograms.
"""
import argparse
import os
import re
import sys
import time

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

    def coverage(self, lx, ly, rotations_rad):
        """Sum over rotations of the maps at camera coords; lx, ly in deg (tangent plane,
        x = east, y = north).  Returns (M, nf)."""
        out = np.zeros((lx.shape[0], self.nf), dtype=np.float32)
        for phi in rotations_rad:
            c, s = np.cos(phi), np.sin(phi)
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


# ----------------------------------------------------------------------------- optimizer
class MapTilingOptimizer:
    def __init__(self, camera, randoms_xyz, tiles_xyz, rotations_deg, max_drift, delta,
                 temp, moves_per_tile, seed=0):
        from scipy.spatial import cKDTree
        self.cam = camera
        self.rot = [np.radians(r) for r in rotations_deg]
        self.rand = np.ascontiguousarray(randoms_xyz, dtype=np.float64)
        self.n_rand = len(self.rand)
        self.tiles = np.ascontiguousarray(tiles_xyz, dtype=np.float64)
        self.tiles0 = self.tiles.copy()
        self.n_tiles = len(self.tiles)
        self.max_drift = np.radians(max_drift)
        self.sigma = np.radians(delta)
        self.temp = temp
        self.moves_per_tile = max(1, int(moves_per_tile))
        self.rng = np.random.default_rng(seed)
        self.cov = np.zeros((self.n_rand, self.cam.nf), dtype=np.float32)
        self.interior = np.ones(self.n_rand, dtype=bool)

        t0 = time.time()
        self.tree = cKDTree(self.rand)
        chord = 2 * np.sin((np.radians(self.cam.max_radius) + self.max_drift) / 2)
        self.cand = self.tree.query_ball_point(self.tiles, chord)
        self.cand = [np.array(c, dtype=np.int64) for c in self.cand]
        print(f'[timing] candidate sets in {time.time() - t0:.1f} s; '
              f'mean {np.mean([len(c) for c in self.cand]):.0f} randoms per tile')

        t0 = time.time()
        for i in range(self.n_tiles):
            c = self.cand[i]
            if len(c):
                self.cov[c] += self._tile_coverage(self.tiles[i], self.rand[c])
        print(f'[timing] initial coverage in {time.time() - t0:.1f} s')
        self._refresh_sums()

    def _tile_coverage(self, t, pts):
        east, north = east_north(t)
        lx, ly = tangent_coords(pts, t, east, north)
        return self.cam.coverage(lx, ly, self.rot)

    def _refresh_sums(self):
        self.sum_c = self.cov.sum(axis=0, dtype=np.float64)
        self.sum_c2 = (self.cov.astype(np.float64) ** 2).sum(axis=0)

    def stats(self):
        mean = self.sum_c / self.n_rand
        var = np.maximum(self.sum_c2 / self.n_rand - mean ** 2, 0)
        return mean, np.sqrt(var)

    def energy(self):
        return float(np.sum(self.sum_c2 - self.sum_c ** 2 / self.n_rand))

    def iterate(self):
        """One sweep over all tiles in random order; returns summary dict."""
        order = self.rng.permutation(self.n_tiles)
        acc = imp = rej = 0
        dists = []
        for i in order:
            c = self.cand[i]
            if len(c) == 0:
                rej += 1
                continue
            pts = self.rand[c]
            cur = self.cov[c].astype(np.float64)
            t = self.tiles[i]
            old = self._tile_coverage(t, pts).astype(np.float64)
            east, north = east_north(t)
            best = None
            for _ in range(self.moves_per_tile):
                dx, dy = self.rng.normal(0, self.sigma, 2)
                p = t + dx * east + dy * north
                p /= np.linalg.norm(p)
                d0 = np.arccos(np.clip(p @ self.tiles0[i], -1, 1))
                if d0 > self.max_drift:
                    p = self.tiles0[i] + (p - self.tiles0[i]) * (self.max_drift / d0)
                    p /= np.linalg.norm(p)
                new = self._tile_coverage(p, pts).astype(np.float64)
                delta = new - old                                   # (C, nf)
                d_sum = delta.sum(axis=0)
                d_sum2 = (2 * cur * delta + delta ** 2).sum(axis=0)
                d_e = float(np.sum(d_sum2 - (2 * self.sum_c * d_sum + d_sum ** 2) / self.n_rand))
                if best is None or d_e < best[0]:
                    best = (d_e, p, delta, d_sum, d_sum2)
            d_e, p, delta, d_sum, d_sum2 = best
            if d_e < 0:
                accept, imp = True, imp + 1
            elif self.temp > 0 and self.rng.random() < np.exp(-d_e / (self.temp * 1e6)):
                accept = True
            else:
                accept = False
            if accept:
                self.cov[c] += delta.astype(np.float32)
                self.sum_c += d_sum; self.sum_c2 += d_sum2
                dists.append(np.degrees(np.arccos(np.clip(p @ t, -1, 1))))
                self.tiles[i] = p
                acc += 1
            else:
                rej += 1
        self._refresh_sums()   # guard against float32 drift
        mean, rms = self.stats()
        return dict(accepted=acc, improved=imp, rejected=rej, mean=mean, rms=rms,
                    move=(np.min(dists), np.median(dists), np.max(dists)) if dists else (0, 0, 0))

    # ----------------------------------------------------------------- output
    def write_tiles(self, path, iteration):
        from astropy.io import fits
        ra, dec = xyz_to_radec(self.tiles); ra0, dec0 = xyz_to_radec(self.tiles0)
        cols = [fits.Column(name='tileid', format='K', array=np.arange(self.n_tiles) + 1),
                fits.Column(name='ra', format='D', array=ra), fits.Column(name='dec', format='D', array=dec),
                fits.Column(name='ra0', format='D', array=ra0), fits.Column(name='dec0', format='D', array=dec0)]
        hdu = fits.BinTableHDU.from_columns(cols, name='TILES')
        hdu.header['ITER'] = iteration
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

    def plot(self, path, ra_center, dec_center, diameter, res_arcmin=0.5):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        ra, dec = xyz_to_radec(self.rand)
        dra = ((ra - ra_center + 180) % 360 - 180) * np.cos(np.radians(dec_center))
        ddec = dec - dec_center
        h = diameter / 2
        sel = (np.abs(dra) < h) & (np.abs(ddec) < h)
        nb = int(diameter * 60 / res_arcmin)
        nf = self.cam.nf
        fig, axes = plt.subplots(2, nf, figsize=(4.6 * nf, 8.5))
        vmax = float(np.percentile(self.cov, 99.5)) * 1.05
        tra, tdec = xyz_to_radec(self.tiles)
        tdra = ((tra - ra_center + 180) % 360 - 180) * np.cos(np.radians(dec_center)); tddec = tdec - dec_center
        tsel = (np.abs(tdra) < h) & (np.abs(tddec) < h)
        for k, f in enumerate(self.cam.filters):
            ax = axes[0, k]
            H, xe, ye = np.histogram2d(dra[sel], ddec[sel], bins=nb, range=[[-h, h], [-h, h]], weights=self.cov[sel, k])
            N, _, _ = np.histogram2d(dra[sel], ddec[sel], bins=nb, range=[[-h, h], [-h, h]])
            im = ax.imshow((H / np.maximum(N, 1)).T, origin='lower', extent=[h, -h, -h, h], vmin=0, vmax=vmax, cmap='viridis')
            ax.plot(tdra[tsel], tddec[tsel], 'w+', ms=4, mew=0.8)
            ax.set_title(f'{f}: coverage around ({ra_center}, {dec_center})'); ax.set_xlabel('dRA cos(dec) [deg]'); ax.set_ylabel('dDec [deg]')
            plt.colorbar(im, ax=ax, fraction=0.046)
            ax = axes[1, k]
            ax.hist(self.cov[:, k], bins=np.arange(0, vmax + 0.25, 0.25), color=f'C{k}')
            mean, rms = self.stats()
            ax.set_title(f'{f}: mean {mean[k]:.2f}, rms {rms[k]:.2f}, rms/mean {rms[k] / mean[k]:.3f}')
            ax.set_xlabel('coverage (sum of throughput)'); ax.set_yscale('log')
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
    p.add_argument('--margin', type=float, default=0.4, help='tiles are placed out to this far beyond the footprint [deg]')
    p.add_argument('--randoms', help='FITS file of randoms (RA, DEC) to reuse')
    p.add_argument('--num-randoms', type=int, default=2_000_000)
    p.add_argument('--tiles', help='initial tiles FITS (ra, dec); a <prefix>tiles_NNNN.fits checkpoint resumes')
    p.add_argument('--ntiles', type=int, help='number of Fibonacci tiles over footprint+margin')
    p.add_argument('--tilearea', type=float, help='Fibonacci tile area [deg^2] (alternative to --ntiles)')
    p.add_argument('--init-offsets', help='start from the footprint pointings with these dither offsets '
                                           '"dx,dy dx,dy ..." [arcmin, +x = +RA]; e.g. Hironao phase 1')
    p.add_argument('--max-drift', type=float, default=1.0, help='max drift from initial position [deg]')
    p.add_argument('--delta', type=float, default=0.1, help='proposal step sigma [deg]')
    p.add_argument('--shrink', type=float, default=1.0, help='multiply delta by this each iteration')
    p.add_argument('--temp', type=float, default=0.0, help='Metropolis temperature')
    p.add_argument('--moves-per-tile', type=int, default=3)
    p.add_argument('--iters', type=int, default=20)
    p.add_argument('--seed', type=int, default=1)
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
    if opts.footprint:
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
    if opts.tiles:
        d = fits.getdata(opts.tiles)
        tiles = radec_to_xyz(d['ra'], d['dec'])
        it = parse_resume_iteration(opts.tiles)
        if it is not None and 'ra0' in d.names:
            start_iter = it
            tiles0 = radec_to_xyz(d['ra0'], d['dec0'])
        else:
            tiles0 = tiles.copy()
        print(f'read {len(tiles)} tiles from {opts.tiles}' + (f', resuming at iteration {it}' if it is not None else ''))
    elif opts.init_offsets:
        off = np.radians(np.array([[float(v) for v in o.split(',')] for o in opts.init_offsets.split()]) / 60.0)
        east, north = east_north(fp.centers)
        tiles = (fp.centers[:, None, :] + off[None, :, 0, None] * east[:, None, :] + off[None, :, 1, None] * north[:, None, :]).reshape(-1, 3)
        tiles /= np.linalg.norm(tiles, axis=1, keepdims=True)
        tiles0 = tiles.copy()
        print(f'{len(tiles)} tiles = {len(fp.centers)} pointings x {len(off)} offsets')
    else:
        if fp is None:
            sys.exit('need --footprint for a Fibonacci initial tiling')
        if opts.ntiles:
            area_ext = fp.area(rng) if opts.margin == 0 else Footprint(fra, fdec, opts.footprint_radius + opts.margin).area(rng)
            tilearea = area_ext / opts.ntiles
        elif opts.tilearea:
            tilearea = opts.tilearea
        else:
            sys.exit('give --ntiles or --tilearea for the Fibonacci initial tiling')
        n_all = int(round(SPHERE_AREA_DEG2 / tilearea))
        _, fra_all, fdec_all = fibonacci_points(n_all)
        # random rotation about the pole so the lattice seam is not aligned with the footprint
        fra_all = (fra_all + rng.uniform(0, 360)) % 360
        xyz = radec_to_xyz(fra_all, fdec_all)
        keep = fp.contains(xyz, extra=opts.margin)
        tiles = xyz[keep]; tiles0 = tiles.copy()
        print(f'Fibonacci lattice: {tilearea:.4f} deg^2 per tile ({n_all} all-sky), {len(tiles)} inside footprint + {opts.margin} deg margin')
    print(f'expected mean coverage per filter ~ {np.round(len(tiles) * len(rot) * cam.area / (fp.area(rng) if fp else 1), 3)} '
          f'(tiles x rotations x area / footprint area; edge tiles cover outside)')

    opt = MapTilingOptimizer(cam, rand, tiles, rot, opts.max_drift, opts.delta, opts.temp, opts.moves_per_tile, seed=opts.seed)
    opt.tiles0 = tiles0
    if fp is not None:
        opt.interior = fp.interior(rand)
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
        print(f'Iter {it:4d} | {s} | sum rms/mean {np.sum(r["rms"] / r["mean"]):.4f} | sigma {np.degrees(opt.sigma):.4f} '
              f'| acc {r["accepted"]} (imp {r["improved"]}) rej {r["rejected"]} | move [{r["move"][0]:.4f}, {r["move"][1]:.4f}, {r["move"][2]:.4f}] '
              f'| {time.time() - t0:.1f} s', flush=True)
        opt.write_tiles(f'{opts.output}tiles_{it:04d}.fits', it)
        if opts.plot:
            opt.plot(f'{opts.output}coverage_{it:04d}.png', opts.ra_center, opts.dec_center, opts.diameter)
        opt.sigma *= opts.shrink
    opt.report('final')


if __name__ == '__main__':
    main()
