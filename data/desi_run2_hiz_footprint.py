"""Sky maps of the DESI Run 2 High-z Dark Energy Survey footprint built from the imaging stripes
(figures/hiz_footprint.png and figures/hiz_stripes_radec.png of the DESI Run 2 document).

Needs Anand Raichoor's footprints-v0.fits (nside 128; $FOOTPRINTS_TABLE, default ~/data/footprints/footprints-v0.fits)
for the E(B-V) background and the Rubin WFD layer.  Run from anywhere: python data/desi_run2_hiz_footprint.py

Layers: SFD E(B-V) background and Rubin WFD from Anand Raichoor's footprints-v0.fits (nside 128);
the imaging stripes from the DESI Run 2 imaging plan (data/desi_run2_stripes.csv), with the Niji
stripes at the bounds adopted for the HSC tiling (data/subaru3/niji_stripes_2p8_hsc.csv).
"""
import csv, os, sys, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from astropy.io import fits
from scipy.spatial import cKDTree

REPO = os.environ.get('SKYTILING_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RA0 = 90.0            # RA at the centre of the map (east to the left, as in the existing figures)
MIN_PIECE_DEG2 = 10.0  # IBIS+POLU pieces of the Niji stripes smaller than this are dropped

def lon(ra):
    """Mollweide longitude [rad] for RA [deg]: RA increases to the left, RA0 at the centre."""
    return -np.radians((np.asarray(ra, float) - RA0 + 180.0) % 360.0 - 180.0)

def read_stripes():
    full = [r for r in csv.DictReader(l for l in open(f'{REPO}/data/desi_run2_stripes.csv') if not l.startswith('#'))]
    niji = {r['name']: r for r in csv.DictReader(l for l in open(f'{REPO}/data/subaru3/niji_stripes_2p8_hsc.csv') if not l.startswith('#'))}
    stripes = []
    for r in full:
        name = r['name']; dc = float(r['dec_center'])
        ext = 'extended' in r['surveys']
        if name in niji and name in ('NGC-5', 'NGC-6', 'SGC-1', 'SGC-2'):
            n = niji[name]
            stripes.append(dict(name=name, ra_min=float(n['ra_min']), ra_max=float(n['ra_max']), dec_min=float(n['dec_min']), dec_max=float(n['dec_max']), kind='niji', ext=False))
        else:
            stripes.append(dict(name=name, ra_min=float(r['ra_min']), ra_max=float(r['ra_max']), dec_min=dc - 1.4, dec_max=dc + 1.4, kind='ibis', ext=ext))
    for name in ('SGC-3', 'SGC-4'):           # the Niji southern extension, in addition to the IBIS+POLU stripes of the same name
        n = niji[name]
        stripes.append(dict(name=name + ' (Niji)', ra_min=float(n['ra_min']), ra_max=float(n['ra_max']), dec_min=float(n['dec_min']), dec_max=float(n['dec_max']), kind='niji', ext=False))
    # the parts of the Niji stripes outside the HSC data, at the DESI RA bounds, are imaged by IBIS+POLU instead
    for r in full:
        name = r['name']
        if name not in ('NGC-5', 'NGC-6', 'SGC-1', 'SGC-2'):
            continue
        n = niji[name]; dc = float(r['dec_center'])
        for side, (lo, hi) in (('W', (float(r['ra_min']), float(n['ra_min']))), ('E', (float(n['ra_max']), float(r['ra_max'])))):
            if (hi - lo) % 360.0 <= 0 or (hi - lo) % 360.0 > 180:
                continue
            piece = dict(name=f'{name} ({side}, IBIS)', ra_min=lo, ra_max=hi, dec_min=dc - 1.4, dec_max=dc + 1.4, kind='ibis', ext=False)
            if rect_area(piece) >= MIN_PIECE_DEG2:
                stripes.append(piece)
    return stripes

def rect_area(s):
    w = (s['ra_max'] - s['ra_min']) % 360.0
    return w * np.degrees(np.sin(np.radians(s['dec_max'])) - np.sin(np.radians(s['dec_min'])))

def rect_poly(s, n=200):
    """Closed polygon (ra, dec) with dense edges, RA unwrapped from ra_min."""
    w = (s['ra_max'] - s['ra_min']) % 360.0
    ra = s['ra_min'] + np.linspace(0, w, n)
    top = (ra, np.full(n, s['dec_max'])); bot = (ra[::-1], np.full(n, s['dec_min']))
    return np.concatenate([top[0], bot[0], [ra[0]]]) % 360.0, np.concatenate([top[1], bot[1], [s['dec_max']]])

def in_rect(ra, dec, s):
    w = (s['ra_max'] - s['ra_min']) % 360.0
    return (((ra - s['ra_min']) % 360.0) <= w) & (dec >= s['dec_min']) & (dec <= s['dec_max'])

def load_table():
    d = fits.open(os.environ.get('FOOTPRINTS_TABLE', os.path.expanduser('~/data/footprints/footprints-v0.fits')))[1].data
    xyz = np.c_[np.cos(np.radians(d['DEC'])) * np.cos(np.radians(d['RA'])), np.cos(np.radians(d['DEC'])) * np.sin(np.radians(d['RA'])), np.sin(np.radians(d['DEC']))]
    return d, cKDTree(xyz)

def grid(res=0.25):
    ra = np.arange(0, 360, res) + res / 2; dec = np.arange(-90, 90, res) + res / 2
    RA, DEC = np.meshgrid(ra, dec)
    return RA, DEC

def lookup(tree, RA, DEC):
    xyz = np.c_[np.cos(np.radians(DEC.ravel())) * np.cos(np.radians(RA.ravel())), np.cos(np.radians(DEC.ravel())) * np.sin(np.radians(RA.ravel())), np.sin(np.radians(DEC.ravel()))]
    _, j = tree.query(xyz); return j.reshape(RA.shape)

def draw_mollweide(ax, d, j, RA, DEC):
    """E(B-V) background, Rubin WFD in blue, graticule."""
    ebv = np.log10(np.clip(d['EBV'][j], 0.005, 2.0))
    L = lon(RA); B = np.radians(DEC)
    # pcolormesh needs monotonic lon: build the mesh in lon order
    order = np.argsort(L[0]); Ls = L[:, order]; Bs = B[:, order]; es = ebv[:, order]
    ax.pcolormesh(Ls, Bs, es, cmap='Greys', vmin=-2.3, vmax=1.4, shading='nearest', rasterized=True, zorder=0)
    rub = d['RUBIN_WFD'][j].astype(float)[:, order]
    ax.contourf(Ls, Bs, rub, levels=[0.5, 1.5], colors=['#8fb8de'], alpha=0.55, zorder=1)
    ax.grid(True, color='0.5', lw=0.3, alpha=0.6)
    ticks = np.array([240, 210, 180, 150, 120, 60, 30, 0, 330, 300])
    ax.set_xticks(lon(ticks)); ax.set_xticklabels([f'{t}°' for t in ticks], fontsize=8)
    ax.set_yticklabels([f'{t:g}°' for t in (-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75)], fontsize=8)
    ax.set_xlabel('R.A. [deg]'); ax.set_ylabel('Dec. [deg]')

def draw_stripes(ax, stripes, RA, DEC, with_outline=True, lw=1.2):
    col = {'ibis': '#f2a65a', 'niji': '#5aa469'}
    for s in stripes:
        ra, dec = rect_poly(s)
        # split at the RA wrap so the polygon does not cross the map edge
        L = lon(ra); B = np.radians(dec)
        kw = dict(facecolor=col[s['kind']], edgecolor='none', alpha=0.85 if not s['ext'] else 0.45, zorder=3)
        if s['ext']:
            kw.update(hatch='////', edgecolor='#b06d1f', lw=0.0)
        brk = np.where(np.abs(np.diff(L)) > np.pi)[0]
        if len(brk) == 0:
            ax.fill(L, B, **kw)
        else:
            # rasterise instead for wrapped rectangles
            m = in_rect(RA, DEC, s).astype(float)
            order = np.argsort(lon(RA)[0]); ax.contourf(lon(RA)[:, order], np.radians(DEC)[:, order], m[:, order], levels=[0.5, 1.5], colors=[col[s['kind']]], alpha=kw['alpha'], zorder=3)
    if with_outline:
        m = np.zeros(RA.shape, bool)
        for s in stripes:
            m |= in_rect(RA, DEC, s)
        order = np.argsort(lon(RA)[0])
        ax.contour(lon(RA)[:, order], np.radians(DEC)[:, order], m[:, order].astype(float), levels=[0.5], colors='r', linewidths=lw, zorder=4)
        return m

def areas(stripes):
    a = {}
    for s in stripes:
        a[s['name']] = rect_area(s)
    ibis_base = sum(rect_area(s) for s in stripes if s['kind'] == 'ibis' and not s['ext'])
    ibis_ext = sum(rect_area(s) for s in stripes if s['kind'] == 'ibis' and s['ext'])
    niji = sum(rect_area(s) for s in stripes if s['kind'] == 'niji')
    niji_own = sum(rect_area(s) for s in stripes if s['kind'] == 'niji' and '(Niji)' not in s['name'])   # Niji stripes not inside an IBIS stripe
    return a, ibis_base, ibis_ext, niji, niji_own

if __name__ == '__main__':
    stripes = read_stripes()
    a, ibis_base, ibis_ext, niji, niji_own = areas(stripes)
    for k, v in a.items(): print(f'{k:14s} {v:7.1f}')
    print(f'IBIS+POLU baseline {ibis_base:.0f}, extended {ibis_ext:.0f}; Niji {niji:.0f} (of which {niji - niji_own:.0f} inside IBIS stripes); high-z baseline {ibis_base + niji_own:.0f}, with extensions {ibis_base + ibis_ext + niji_own:.0f}')
    d, tree = load_table()
    RA, DEC = grid(0.25); j = lookup(tree, RA, DEC)
    rub_area = d['RUBIN_WFD'].sum() * 4 * np.pi * (180 / np.pi) ** 2 / len(d)
    fig = plt.figure(figsize=(10, 6.2)); ax = fig.add_subplot(111, projection='mollweide')
    draw_mollweide(ax, d, j, RA, DEC)
    m = draw_stripes(ax, stripes, RA, DEC)
    cell = np.radians(0.25) ** 2 * np.cos(np.radians(DEC)) * (180 / np.pi) ** 2
    print('union of stripes (raster) %.0f deg2' % (cell * m).sum())
    handles = [Patch(facecolor='#8fb8de', alpha=0.55, label=f'Rubin/WFD ({rub_area/1000:.0f}k deg$^2$)'),
               Patch(facecolor='#f2a65a', alpha=0.85, label=f'IBIS+POLU stripes ({ibis_base:,.0f} deg$^2$)'),
               Patch(facecolor='#f2a65a', alpha=0.45, hatch='////', edgecolor='#b06d1f', label=f'possible extensions ({ibis_ext:,.0f} deg$^2$)'),
               Patch(facecolor='#5aa469', alpha=0.85, label=f'HSC-Niji+VST stripes ({niji:,.0f} deg$^2$)'),
               Line2D([], [], color='r', lw=1.2, label=f'Imaging footprint = High-z Dark Energy Survey ({ibis_base + niji_own:,.0f} deg$^2$; {ibis_base + ibis_ext + niji_own:,.0f} with extensions)')]
    ax.legend(handles=handles, loc='lower center', fontsize=8, framealpha=0.95, ncol=2, bbox_to_anchor=(0.5, -0.22))
    plt.tight_layout(); plt.savefig('hiz_footprint.png', dpi=150); plt.savefig('hiz_footprint.pdf')
    # zoom panels in RA/Dec for the stripes themselves
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2))
    for ax, (lab, lo, hi) in zip(axes, (('NGC', 115, 260), ('SGC', 295, 440))):
        for s in stripes:
            if not s['name'].startswith(lab): continue
            ra, dec = rect_poly(s, 2); x0 = s['ra_min'] if lab == 'NGC' else (s['ra_min'] - 295) % 360 + 295
            w = (s['ra_max'] - s['ra_min']) % 360.0
            kw = dict(facecolor='#f2a65a' if s['kind'] == 'ibis' else '#5aa469', alpha=0.85 if not s['ext'] else 0.45, edgecolor='k', lw=0.4)
            if s['ext']: kw.update(hatch='////', edgecolor='#b06d1f')
            ax.add_patch(plt.Rectangle((x0, s['dec_min']), w, s['dec_max'] - s['dec_min'], **kw))
            if 'IBIS)' not in s['name']: ax.text(x0 + w + 0.6, 0.5 * (s['dec_min'] + s['dec_max']), s['name'].replace(' (Niji)', ''), fontsize=7, va='center')
        ax.set_xlim(hi, lo); ax.set_ylim(-23, 18); ax.grid(alpha=0.3)
        ticks = np.arange(lo, hi + 1, 15); ax.set_xticks(ticks); ax.set_xticklabels([f'{t % 360:g}' for t in ticks])
        ax.set_xlabel('R.A. [deg]'); ax.set_ylabel('Dec. [deg]'); ax.set_title(f'{lab}: DESI Run 2 imaging stripes (orange: IBIS+POLU, hatched: possible extensions; green: HSC-Niji+VST)', fontsize=9)
    plt.tight_layout(); plt.savefig('stripes_radec.png', dpi=150); plt.savefig('stripes_radec.pdf')
    print('done')
