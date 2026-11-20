#!/usr/bin/env python3
"""Produce separate figures using only measured CSV, with explicit Euclidean gauge."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'results/phase_fixed_plane_profile.csv').open() as f: rows=list(csv.DictReader(f))
std=[r for r in rows if r['tolerance']=='standard'];tight=[r for r in rows if r['tolerance']=='tight']
x=np.asarray([float(t['xi']) for t in std]);xt=np.asarray([float(t['xi']) for t in tight])
fig,ax=plt.subplots(figsize=(10,5))
for label,k in [('total','grassmann_total_speed_euclidean_gauge'),('transport','grassmann_transport_speed_euclidean_gauge'),('source birth','grassmann_birth_speed_euclidean_gauge')]:
 ax.plot(x,[float(t[k]) for t in std],label=label)
ax.set(xlabel='xi',ylabel='Grassmann speed (10D Euclidean gauge)',title='v18.11 — Material 2-plane rotation: homogeneous transport vs source birth')
ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'results/grassmann_speed_transport_birth.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5))
for label,table,XX,alpha in [('standard',std,x,1),('tight',tight,xt,.5)]:
 ax.plot(XX,[float(t['plane_next_step_max_principal_angle_rad_euclidean_gauge']) for t in table],label=label,alpha=alpha)
ax.set(xlabel='xi',ylabel='largest successive-plane principal angle (rad)',title='v18.11 — Inter-sample plane motion; 10D Euclidean gauge')
ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'results/plane_principal_increment.png',dpi=170);plt.close(fig)
print('Wrote two v18.11 diagnostic figures.')
