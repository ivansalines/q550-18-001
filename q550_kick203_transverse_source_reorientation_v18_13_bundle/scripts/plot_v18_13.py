#!/usr/bin/env python3
"""Plots produced from current v18.13 CSV, no external figure provenance."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results'

def rows(path):
    with (R/path).open() as f:return list(csv.DictReader(f))

s=[q for q in rows('source_components_profile.csv') if q['tolerance']=='standard']
t=[q for q in rows('source_components_profile.csv') if q['tolerance']=='tight']
e=[q for q in rows('event_turnings.csv') if q['tolerance']=='standard']
xx=[float(q['xi']) for q in s];tt=[float(q['xi']) for q in t]
fig,ax=plt.subplots(figsize=(10,4.8))
ax.plot(xx,[float(q['radial_rate']) for q in s],label='radial source rate — standard',linewidth=1.6)
ax.plot(tt,[float(q['radial_rate']) for q in t],label='radial source rate — tight',alpha=.5,linewidth=.8,linestyle=':')
ax.axhline(0,color='gray',lw=.6)
for z in e:ax.axvline(float(z['xi_pchip']),alpha=.27,lw=1,linestyle='--')
ax.set(xlabel='xi',ylabel='d rho / d xi',title='Radial source transport: four signed turnings')
ax.legend(loc='upper right');fig.tight_layout();fig.savefig(R/'radial_source_turnings.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,4.8))
ax.plot(xx,[float(q['turn_speed']) for q in s],label='residual directional speed — standard',linewidth=1.7)
ax.plot(tt,[float(q['turn_speed']) for q in t],label='tight',lw=.8,linestyle=':')
for z in e:ax.axvline(float(z['xi_pchip']),alpha=.25,linestyle='--',lw=1)
ax.set(xlabel='xi',ylabel='||dv/dxi||  (initial Euclidean gauge)',title='Residual plane continues to turn through radial reversals')
ax.legend();fig.tight_layout();fig.savefig(R/'residual_direction_speed_events.png',dpi=160);plt.close(fig)
w=[q for q in rows('event_orientation_windows.csv') if q['tolerance']=='standard']
fig,ax=plt.subplots(figsize=(9,4.7))
for i in range(1,5):
  q=[a for a in w if int(a['event'])==i]
  ax.plot([float(a['halfwidth_xi']) for a in q],[float(a['cosine']) for a in q],'-o',label=f'event {i} @ xi={float(q[0]["xi_pchip"]):.3f}')
ax.axhline(0,lw=.7,color='gray')
ax.set(xlabel='half-window in xi',ylabel='parallel-transported tangent cosine',ylim=(-1.08,1.08),
       title='Directional coherence across radial turning: local vs wide windows')
ax.legend();fig.tight_layout();fig.savefig(R/'event_orientation_coherence.png',dpi=160);plt.close(fig)
print('plots complete')
