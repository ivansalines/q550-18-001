"""Plots genuine frozen rank spectrum and, when available, server nonlinear validation."""
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results'
with open(out/'jet_rank_hierarchy.csv') as f:rows=list(csv.DictReader(f))
fig,axs=plt.subplots(1,3,figsize=(13,3.9))
for ax,name in zip(axs,['material3','raw4','raw6']):
 s=[r for r in rows if r['xi']=='4.0' and r['observer']==name and r['tau']=='0.1']
 ax.plot([int(r['jet_order']) for r in s],[int(r['rank']) for r in s],'-o',linewidth=2)
 ax.set(title=name,xlabel='Highest jet order',ylabel='Observable rank',ylim=(0,11),xticks=[0,1,2,3]);ax.grid(alpha=.3)
fig.suptitle('Q550 v18.20 — local 10D observability hierarchy at xi=4')
fig.tight_layout();fig.savefig(out/'jet_rank_hierarchy.png',dpi=160);plt.close(fig)
p=out/'nonlinear_jet_validation.csv'
if p.exists():
 with p.open() as f:rr=list(csv.DictReader(f))
 fig,ax=plt.subplots(figsize=(8,4.5))
 for label in sorted(set(r['label'] for r in rr)):
  ss=[r for r in rr if r['label']==label and r['epsilon']=='3e-06' and float(r['delta_xi'])>0]
  if ss:ax.plot([float(r['delta_xi']) for r in ss],[float(r['scaled_reduced_separation']) for r in ss],'-o',label=label)
 ax.set(xlabel='xi - 4',ylabel='Reduced gap / 2 epsilon',title='Independent nonlinear paired futures (Q550 original ODE)')
 ax.legend(fontsize=7);ax.grid(alpha=.3);fig.tight_layout();fig.savefig(out/'nonlinear_jet_growth.png',dpi=160);plt.close(fig)
