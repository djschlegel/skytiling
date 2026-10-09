import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys; sys.path.insert(0,'/home/claude/hsc/stripes_eval')
a=np.load('evalA.npy',allow_pickle=True).item(); b=np.load('evalB.npy',allow_pickle=True).item(); a.update(b)
runs=[('ours_p2','skytiling stripes, rms/mean$^2$ (10,880 exp.)','C3'),('hironao_rs4m5_l16','Hironao rs4m5_l16, var + 16 x holes (12,832 exp.)','C0'),('hironao_init','Hironao start: 802 pointings x 4+2 square x 4 rot','0.5')]
fig,axes=plt.subplots(2,2,figsize=(14,8))
for j,lab in enumerate(('NGC','SGC')):
    for k,(name,title,col) in enumerate(runs):
        x,p=a[name]['prof'][lab]
        axes[0,j].plot(x,p[:,0],color=col,lw=1.2 if k<2 else 0.8,ls='-' if k<2 else '--',label=title)
        axes[1,j].plot(x,p[:,4],color=col,lw=1.2 if k<2 else 0.8,ls='-' if k<2 else '--')
    edges={'NGC':(-2.12,0.68,3.48),'SGC':(-6.51,-3.71,-0.91,1.89,4.69)}[lab]
    for ax in axes[:,j]:
        for e in edges: ax.axvline(e,color='r',lw=0.6,alpha=0.6)
        ax.grid(alpha=0.3); ax.set_xlabel('Dec [deg]')
    axes[0,j].set_title(f'{lab} band (RA {"140-215" if lab=="NGC" else "340-28"}): mean MBQ1-413 depth vs Dec (red: stripe edges)',fontsize=10)
    axes[0,j].set_ylabel('mean depth (sum of THRU), 413'); axes[1,j].set_ylabel('mean n_exp (GOOD), 413')
axes[0,0].legend(fontsize=8,loc='lower center')
plt.tight_layout(); plt.savefig('compare_dec_profiles.png',dpi=110)
# histograms of 413 depth, stripe interior
fig,ax=plt.subplots(figsize=(7,4.5)); bins=np.arange(0,8.01,0.1); xc=0.5*(bins[1:]+bins[:-1])
for name,title,col in runs:
    h=a[name]['hist']['stripes_int']; ax.step(xc,h/h.sum()/0.1,where='mid',color=col,label=title,lw=1.2)
ax.set_xlabel('MBQ1-413 depth (sum of THRU) in the stripe interior'); ax.set_ylabel('fraction per unit depth'); ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_xlim(0,7)
plt.tight_layout(); plt.savefig('compare_hist413.png',dpi=110); print('ok')
