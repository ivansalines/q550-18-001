#!/usr/bin/env python3
"""Diagnostic visualizations using frozen tabular results."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];RESULT=ROOT/'results'
with (RESULT/'source_kernel_decomposition.csv').open() as f:
    rows=[r for r in csv.DictReader(f) if r['tolerance']=='standard' and r['metric']=='euclidean']
xs=np.array(sorted(set(float(r['xi']) for r in rows)))
ls=np.array(sorted(set(float(r['lag']) for r in rows)))
fields=['exact','instantaneous','amplitude_gated_unit','covariance_residual']
arrays={k:np.array([float(r[k]) for r in rows]).reshape(len(xs),len(ls)) for k in fields}
fig,axs=plt.subplots(2,2,figsize=(12,7),sharex=True,sharey=True,constrained_layout=True)
for ax,key,name in zip(axs.flat,fields,['Exact source-integrated K','Endpoint instantaneous source','Amplitude-gated unit pairs','Orientation-amplitude covariance']):
    im=ax.pcolormesh(xs,ls,arrays[key].T,shading='auto',cmap='coolwarm',vmin=-1 if key!='covariance_residual' else None,vmax=1 if key!='covariance_residual' else None)
    fig.colorbar(im,ax=ax,label='Kernel / diagnostic')
    ax.set_title(name);ax.set_xlabel('xi');ax.set_ylabel('Lag (xi)')
fig.suptitle('v18.18 — source-only closure vs incomplete recurrence observers')
fig.savefig(RESULT/'source_recurrence_residual_maps.png',dpi=140);plt.close(fig)
with (RESULT/'time_residual_and_ridge.csv').open() as f:
    r=[a for a in csv.DictReader(f) if a['tolerance']=='standard' and a['metric']=='euclidean']
x=np.array([float(v['xi']) for v in r]); fig,ax=plt.subplots(figsize=(10,4))
for field,label in [('lag_exact','Integrated kernel'),('lag_instantaneous','Instantaneous endpoint'),('lag_unit_only','Unit direction pairs')]:
    ax.plot(x,[float(v[field]) for v in r],label=label)
ax.axvline(6.,ls='--',label='Holdout split')
ax.set(xlabel='xi',ylabel='Grid minimum lag (xi)',title='Observed moving ridge and source-only simplifications');ax.legend()
fig.tight_layout();fig.savefig(RESULT/'source_ridge_baseline_comparison.png',dpi=140);plt.close(fig)
