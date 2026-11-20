#!/usr/bin/env python3
"""Reproducible visualizations of the signed source-lag kernel."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results'
with (R/'kernel_xi_lag_profile.csv').open() as f:rows=list(csv.DictReader(f))
profile=[v for v in rows if v['tolerance']=='standard' and v['history']=='fixed_3xi' and v['kernel']]
xs=sorted(set(float(r['xi']) for r in profile));lags=sorted(set(float(r['lag']) for r in profile))
full=np.full((len(lags),len(xs)),np.nan)
xi={v:i for i,v in enumerate(xs)};li={v:i for i,v in enumerate(lags)}
for r in profile:full[li[float(r['lag'])],xi[float(r['xi'])]]=float(r['kernel'])
fig,ax=plt.subplots(figsize=(10,5.2));m=ax.imshow(full,origin='lower',aspect='auto',extent=(xs[0],xs[-1],lags[0],lags[-1]),vmin=-1,vmax=1)
fig.colorbar(m,ax=ax,label='Signed normalized source-lag kernel')
ax.set_xlabel('Observation xi');ax.set_ylabel('Source-age lag (xi)');ax.set_title('Q550 v18.16 — fixed 3 xi lookback, post-anchor source')
ax.axvline(3.9496098,linestyle='--',linewidth=1,alpha=.8);ax.axvline(6.9582933,linestyle='--',linewidth=1,alpha=.8)
fig.tight_layout();fig.savefig(R/'source_lag_kernel_heatmap.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,4.8))
for lag in [1.25,1.5,1.75]:
 rr=sorted((float(r['xi']),float(r['kernel'])) for r in profile if abs(float(r['lag'])-lag)<1e-9)
 ax.plot([v[0] for v in rr],[v[1] for v in rr],label=f'lag={lag:g} xi')
ax.axhline(0,linestyle=':',linewidth=1)
for z in [3.9496098,6.9582933]:ax.axvline(z,linestyle='--',linewidth=.8)
ax.set_xlabel('Observation xi');ax.set_ylabel('Normalized signed kernel')
ax.set_title('Echo persistence: compared lags in a matched recent window');ax.legend();ax.grid(alpha=.18)
fig.tight_layout();fig.savefig(R/'echo_persistence_comparison.png',dpi=160);plt.close(fig)
