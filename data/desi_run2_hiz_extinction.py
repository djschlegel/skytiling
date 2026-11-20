"""E(B-V) statistics and map inside the DESI Run 2 high-z stripes, from the full-resolution SFD maps
(sfd.py; $SFD_MAPS, default ~/dustpub/maps).  Run from anywhere: python data/desi_run2_hiz_extinction.py"""
import numpy as np, os, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
REPO = os.environ.get('SKYTILING_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(REPO, 'data', 'desi_run2_hiz_footprint.py')).read().split('if __name__')[0])
from sfd import sfd_ebv
stripes=read_stripes()
res=0.05
ra=np.arange(0,360,res)+res/2; dec=np.arange(-23,17,res)+res/2; RA,DEC=np.meshgrid(ra,dec)
union=np.zeros(RA.shape,bool)
for s in stripes: union|=in_rect(RA,DEC,s)
base=np.zeros(RA.shape,bool)            # the baseline footprint: without bonus regions (SGC-1)
for s in stripes:
    if not s['ext']: base|=in_rect(RA,DEC,s)
bonus=np.zeros(RA.shape,bool)           # bonus regions outside the baseline (SGC-1, HSC-Niji+VST only)
for s in stripes:
    if s['ext']: bonus|=in_rect(RA,DEC,s)
ebv_all=sfd_ebv(RA.ravel(),DEC.ravel()).reshape(RA.shape)      # the whole panel area
ebv=np.where(union, ebv_all, np.nan)
cell=np.radians(res)**2*np.cos(np.radians(DEC))*(180/np.pi)**2
print('%-16s %7s %6s %6s %6s | area > 0.05 / 0.10 / 0.15 / 0.20 [deg2]'%('stripe','area','median','mean','p90'))
tot=np.zeros(5)
for s in stripes:
    m=in_rect(RA,DEC,s); a=(cell*m).sum(); e=ebv[m]; w=cell[m]
    hi=[(w*(e>t)).sum() for t in (0.05,0.10,0.15,0.20)]
    print('%-16s %7.1f %6.3f %6.3f %6.3f | %6.1f %6.1f %6.1f %6.1f'%(s['name'],a,np.median(e),np.average(e,weights=w),np.percentile(e,90),*hi))
    if not s['ext'] and '(Niji)' not in s['name']: tot+=np.array([a]+hi)
print('baseline total    %7.1f                       | %6.1f %6.1f %6.1f %6.1f'%tuple(tot))
e=ebv[base]; w=cell[base]
for t in (0.05,0.08,0.10,0.15,0.20): print('whole footprint: area with E(B-V) > %.2f: %6.1f deg2 (%.1f%%)'%(t,(w*(e>t)).sum(),100*(w*(e>t)).sum()/w.sum()))
# galactic latitude of the grid, for the |b| = 20 deg lines
from astropy.coordinates import SkyCoord
import astropy.units as u
gb = SkyCoord(ra=RA.ravel()*u.deg, dec=DEC.ravel()*u.deg).galactic.b.deg.reshape(RA.shape)
np.save('ebv_grid.npy', dict(RA=RA,DEC=DEC,ebv=ebv,union=union), allow_pickle=True)
fig, axes = plt.subplots(2, 1, figsize=(14, 7.2))
for ax,(lab,lo,hi) in zip(axes,(('NGC',115,260),('SGC',295,440))):
    cols=np.where((ra>=lo)&(ra<=hi))[0] if lab=='NGC' else np.r_[np.where(ra>=295)[0], np.where(ra<=80)[0]]
    sub=ebv_all[:,cols]; usub=base[:,cols].astype(float); bsub=gb[:,cols]
    im=ax.imshow(sub, origin='lower', extent=[lo,hi,-23,17], aspect='auto', cmap='magma_r', vmin=0, vmax=0.2, interpolation='nearest')
    x=np.linspace(lo,hi,sub.shape[1]); y=np.linspace(-23,17,sub.shape[0])
    ax.contour(x,y,np.where(usub>0.5,sub,0),levels=[0.1],colors='c',linewidths=0.5)
    ax.contour(x,y,usub,levels=[0.5],colors='k',linewidths=2.0)
    if bonus[:,cols].any(): ax.contour(x,y,bonus[:,cols].astype(float),levels=[0.5],colors='k',linewidths=1.2,linestyles='--')   # bonus region (SGC-1)
    ax.contour(x,y,bsub,levels=[-20,20],colors='w',linewidths=1.2,linestyles='--')
    ax.contour(x,y,bsub,levels=[-20,20],colors='k',linewidths=0.4,linestyles='--')
    ax.set_xlim(hi,lo); ax.set_ylim(-23,17); ax.grid(alpha=0.3)
    ticks=np.arange(lo,hi+1,15); ax.set_xticks(ticks); ax.set_xticklabels([f'{t%360:g}' for t in ticks])
    ax.set_xlabel('R.A. [deg]'); ax.set_ylabel('Dec. [deg]')
    ax.set_title(f'{lab}: SFD E(B-V); thick line: High-z Dark Energy Survey footprint (dashed: SGC-1 bonus region); cyan: E(B-V) = 0.1 inside it; dashed: Galactic latitude b = ±20°', fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label='E(B-V)')
plt.tight_layout(); plt.savefig('hiz_stripes_ebv.png', dpi=130); plt.savefig('hiz_stripes_ebv.pdf', dpi=150); print('ok')
