import numpy as np, re, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = np.load('all_eval.npy', allow_pickle=True).item()
names = [('start', 'SSP pointings + dithers (start)', '0.4', ':'), ('ssp_p1', 'SSP start, rms/mean (p=1)', 'C0', '-'), ('final_p2', 'SSP start, rms/mean$^2$ (p=2, adopted)', 'C3', '-'), ('lattice_p1', 'lattice start, p=1', 'C2', '--')]
filters = ['413', '439', '465', '490']
# Fig 1: Dec profiles across the two bands
fig, axes = plt.subplots(2, 2, figsize=(14, 8))
for j, (lab, lo, hi) in enumerate((('NGC band (RA 140-215)', -2.12, 3.48), ('SGC band (RA 340-28)', -0.91, 4.69))):
    for i, k in enumerate((0, None)):
        ax = axes[i, j]
        for n, L, c, ls in names:
            x, p = d[n]['prof'][lab[:3]]
            y = p[:, 0] if k == 0 else p.mean(1)
            ax.plot(x, y, ls, color=c, lw=1.6, label=L)
        for e in (lo, hi): ax.axvline(e, color='k', lw=0.8)
        ax.axvspan(lo, hi, color='0.93', zorder=0)
        ax.set_xlabel('Dec [deg]'); ax.set_ylabel('mean coverage ' + ('MBQ1-413' if k == 0 else 'mean of 4 sub-filters')); ax.grid(alpha=0.3)
        ax.set_title(f'{lab}: coverage profile across the stripe band', fontsize=10)
        if i == 0 and j == 0: ax.legend(fontsize=8, loc='lower center')
plt.tight_layout(); plt.savefig('hsc_stripes_edge_profiles.png', dpi=110); plt.close(fig)
# Fig 2: convergence + histogram
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axes[0]
for n, L, c, ls in names[1:]:
    it, s = [], []
    for line in open(f'{n}.out'):
        m = re.match(r'Iter\s+(\d+) .*sum rms/mean ([\d.]+)', line)
        if m: it.append(int(m.group(1))); s.append(float(m.group(2)))
    ax.plot(it, s, ls, color=c, label=L)
ax.set_xlabel('iteration'); ax.set_ylabel('sum over filters of rms/mean (training randoms, inside footprint)'); ax.grid(alpha=0.3); ax.legend(fontsize=8); ax.set_title('convergence (training randoms; adopted run: 200 iterations on 4M randoms, others 100 on 2M)', fontsize=10)
ax = axes[1]; bins = np.arange(0, 7.01, 0.25)
for n, L, c, ls in names:
    cv = d[n]['cov_in']; ax.hist(cv, bins=bins, histtype='step', color=c, ls=ls, lw=1.5, label=f'{L}: {cv.mean():.2f} ± {cv.std():.2f}')
ax.set_xlabel('MBQ1-413 coverage inside the stripes'); ax.set_ylabel('randoms'); ax.legend(fontsize=8); ax.set_title('coverage distribution, worst sub-filter', fontsize=10)
plt.tight_layout(); plt.savefig('hsc_stripes_stats.png', dpi=110); plt.close(fig)
# table
print('%-12s %-22s %-28s %-28s %-22s' % ('run', 'spill outside [%]', 'inside rms/mean', 'interior rms/mean', 'interior <1.5 [%]'))
for n, L, c, ls in names + [('ssp_p2', '', '', ''), ('cloud_ssp_p2', '', '', '')]:
    if n not in d: continue
    r = d[n]; m, s, h15, h05 = r['inside']; mi, si, hi15, hi05 = r['interior']
    print('%-12s %-22s %-28s %-28s %-22s mean inside %s' % (n, np.round(r['spill'], 1), np.round(s/m, 3), np.round(si/mi, 3), np.round(hi15, 1), np.round(m, 2)))
