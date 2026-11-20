#!/usr/bin/env python3
"""Q550 v18.17 graphics from auditable CSV outputs."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results'
with (R/'moving_ridge_surface.csv').open() as f: surface=list(csv.DictReader(f))
with (R/'ridge_tracking.csv').open() as f: rr=list(csv.DictReader(f))

def view(row,kind,metric='euclidean'):
    return row['tolerance']=='standard' and row['metric']==metric and row['method']==kind

src=[r for r in surface if view(r,'matched_1xi')]
xs=sorted(set(float(s['xi']) for s in src));lags=sorted(set(float(s['lag']) for s in src))
arr=np.full((len(lags),len(xs)),np.nan)
xiidx={x:i for i,x in enumerate(xs)};li={l:j for j,l in enumerate(lags)}
for s in src:arr[li[float(s['lag'])],xiidx[float(s['xi'])]]=float(s['kernel'])
fig,ax=plt.subplots(figsize=(10.5,5.4))
mesh=ax.pcolormesh(xs,lags,arr,shading='nearest',vmin=-1,vmax=1,cmap='coolwarm')
fig.colorbar(mesh,ax=ax,label='Signed recurrence K in fixed paired-age support')
for metric,style in [('euclidean','-'),('strong_diag','--')]:
    curve=sorted((float(s['xi']),float(s['lag_quad'])) for s in rr if view(s,'matched_1xi',metric))
    ax.plot([x for x,y in curve],[y for x,y in curve],linestyle=style,linewidth=2,label=f'{metric} ridge')
ax.set(xlabel='Observation ξ',ylabel='Source-age separation L (ξ)',title='Q550 v18.17 — mobile opposition ridge with fixed 1 ξ paired support')
ax.legend(loc='lower right');fig.tight_layout();fig.savefig(R/'moving_ridge_heatmap.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(10.5,4.9))
for kind in ['matched_1xi','instantaneous_pair','variable_support_3xi','recent_half_xi','older_half_xi']:
    q=sorted((float(s['xi']),float(s['lag_quad'])) for s in rr if view(s,kind))
    ax.plot([v[0] for v in q],[v[1] for v in q],label=kind,linewidth=1.7)
ax.set(xlabel='Observation ξ',ylabel='Refined most-negative lag (ξ)',title='Fixed-support ridge vs endpoint source recurrence and variable overlap')
ax.grid(alpha=.2);ax.legend(loc='best',fontsize=8);fig.tight_layout();fig.savefig(R/'ridge_methods_comparison.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(10,4.8))
for lag in [1.25,1.375,1.5,1.75]:
    j=int(np.argmin(abs(np.array(lags)-lag)))
    ax.plot(xs,arr[j,:],label=f'L={lags[j]:.3f} ξ')
ax.set(xlabel='Observation ξ',ylabel='Signed fixed-support K',title='Persistence of selected source recurrences after removing lag-dependent support')
ax.axhline(0,linestyle=':',linewidth=1);ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(R/'fixed_support_echo_persistence.png',dpi=170);plt.close(fig)
