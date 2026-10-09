"""HSC-Niji rectangles against the Rubin v5.3.5 visit-count maps (rubin_visit_maps.py): fraction of each rectangle
above visit thresholds, Dec profiles at the SGC WFD edge, and the figure data/subaru3/niji_sgc_vs_rubin_v535.png."""
import numpy as np, csv, os, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
out = np.load('rubin_visits_grid.npy', allow_pickle=True).item()
niji = [r for r in csv.DictReader(l for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'subaru3', 'niji_stripes_2p8_hsc.csv')) if not l.startswith('#'))]
def in_rect(RA, DEC, r):
    w = (float(r['ra_max']) - float(r['ra_min'])) % 360
    return (((RA - float(r['ra_min'])) % 360) <= w) & (DEC >= float(r['dec_min'])) & (DEC <= float(r['dec_max']))
print('%-7s %7s | fraction of area with 10-yr visits >= threshold:  u>=30  u>=40  g>=30  g>=45  r>=100  all ugriz >= 30' % ('rect', 'area'))
for r in niji:
    lab = 'NGC' if r['name'].startswith('NGC') else 'SGC'
    g = out[lab]; RA, DEC = g['RA'], g['DEC']; m = in_rect(RA, DEC, r)
    cell = np.radians(0.1)**2 * np.cos(np.radians(DEC)) * (180/np.pi)**2
    a = (cell*m).sum()
    def frac(cond): return 100*(cell*(m & cond)).sum()/a
    allok = (g['u'] >= 30) & (g['g'] >= 30) & (g['r'] >= 30) & (g['i'] >= 30) & (g['z'] >= 30)
    print('%-7s %7.1f | %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%%' % (r['name'], a, frac(g['u'] >= 30), frac(g['u'] >= 40), frac(g['g'] >= 30), frac(g['g'] >= 45), frac(g['r'] >= 100), frac(allok)))
# Dec profiles of u, g, r counts at a few RAs in the SGC
g = out['SGC']; RA, DEC = g['RA'], g['DEC']
print('\n10-yr visits vs Dec in the SGC (u / g / r), at RA 340, 0, 20, 35:')
for dec in np.arange(1.0, 5.6, 0.5):
    row = np.argmin(abs(DEC[:, 0] - dec)); s = '  Dec %+.1f:' % dec
    for ra0 in (340, 0, 20, 35):
        col = np.argmin(abs(((RA[0] - ra0 + 180) % 360) - 180)); s += '   %3d/%3d/%3d' % (g['u'][row, col], g['g'][row, col], g['r'][row, col])
    print(s)
# figure
fig, axes = plt.subplots(2, 1, figsize=(13, 8))
for ax, b, vmax in zip(axes, ('u', 'r'), (70, 200)):
    img = g[b].astype(float); x = ((RA[0] - 320) % 360) + 320
    im = ax.imshow(img, origin='lower', extent=[320, 410, -10, 10], aspect='auto', cmap='viridis', vmin=0, vmax=vmax, interpolation='nearest')
    ax.contour(x, DEC[:, 0], img, levels=[30], colors='w', linewidths=1.0)
    for r in niji:
        if not r['name'].startswith('SGC'): continue
        x0 = float(r['ra_min']); x0 = x0 if x0 > 180 else x0 + 360; w = (float(r['ra_max']) - float(r['ra_min'])) % 360
        ax.add_patch(plt.Rectangle((x0, float(r['dec_min'])), w, float(r['dec_max'])-float(r['dec_min']), fill=False, ec='r', lw=1.5))
        ax.text(x0+0.5, float(r['dec_min'])+0.25, r['name'], color='r', fontsize=9)
    ax.set_xlim(410, 320); ticks = np.arange(320, 411, 15); ax.set_xticks(ticks); ax.set_xticklabels([f'{t%360:g}' for t in ticks])
    ax.set_xlabel('RA [deg]'); ax.set_ylabel('Dec [deg]'); ax.grid(alpha=0.3)
    ax.set_title(f'Rubin baseline v5.3.5, 10-year {b}-band visits (DD/ToO/twilight excluded; white contour: 30 visits); HSC-Niji rectangles in red', fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label=f'{b} visits')
plt.tight_layout(); plt.savefig('niji_sgc_vs_rubin_v535.png', dpi=120); print('ok')
