#!/usr/bin/env python3
from pathlib import Path
import csv
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]; p=root/'results/residual_plane_profile.csv'
rows=list(csv.DictReader(p.open()))
plt.figure(figsize=(10,5))
for tol in ('standard','tight'):
 d=[r for r in rows if r['tolerance']==tol];x=[float(z['xi']) for z in d];y=[float(z['residual_plane_turn_speed_gauge']) for z in d]
 plt.plot(x,y,label=tol,linewidth=1.8 if tol=='standard' else 1,alpha=.85 if tol=='standard' else .7)
for x in (4.135199,4.997337,7.174436):plt.axvline(x,linestyle=':',alpha=.45)
plt.xlabel('xi');plt.ylabel('Residual plane turn rate / xi (initial Euclidean 10D gauge)')
plt.title('Q550 kick203: movement remaining after inverse homogeneous transport')
plt.grid(alpha=.2);plt.legend();plt.tight_layout()
plt.savefig(root/'results/residual_plane_turn_rate.png',dpi=180);plt.close()
plt.figure(figsize=(10,5))
for tol in ('standard','tight'):
 d=[r for r in rows if r['tolerance']==tol];x=[float(z['xi']) for z in d];y=[float(z['cumulative_residual_turn_length_to_here_gauge']) for z in d]
 plt.plot(x,y,label=tol)
plt.xlabel('xi');plt.ylabel('Cumulative turn length (radians, Euclidean initial gauge)')
plt.title('Accumulated path of co-moving plane vs net endpoint displacement')
plt.grid(alpha=.2);plt.legend();plt.tight_layout();plt.savefig(root/'results/residual_plane_cumulative_path.png',dpi=180);plt.close()
