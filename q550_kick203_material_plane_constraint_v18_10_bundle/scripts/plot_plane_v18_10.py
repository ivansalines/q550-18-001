#!/usr/bin/env python3
"""Reproducible static figures: all normalized magnitudes are gauge-dependent."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'results/material_plane_profile.csv').open()))
summary=json.loads((ROOT/'results/summary.json').read_text())
for metric,filename,title in [
 ('source_frame_sine_metric_choice','source_frame_plane_sine.png','Source-frame plane separation (initial Euclidean gauge)'),
 ('source_frame_3volume_gauge_ratio','source_escape_3volume.png','Incoming-source 3-volume (initial Euclidean gauge)'),
 ('lab_sine_metric_choice','laboratory_plane_sine.png','Laboratory plane separation (11D Euclidean gauge)')]:
 fig,ax=plt.subplots(figsize=(9,4.8))
 for tol in ['standard','tight']:
  x=np.array([float(r['xi']) for r in rows if r['tolerance']==tol]);
  y=np.array([float(r[metric]) for r in rows if r['tolerance']==tol]);
  ax.plot(x,y,label=tol,linewidth=1.55 if tol=='standard' else .85,alpha=.85)
 if metric=='lab_sine_metric_choice':
  for root in summary['by_tolerance'][0]['observed_euclidean_dot_crossings_xi']:
   ax.axvline(root,linestyle=':',alpha=.5)
 ax.set(xlabel='xi after kick 203',ylabel='normalized magnitude (not an intrinsic angle)',title=title)
 ax.set_xlim(1,8);ax.grid(alpha=.22);ax.legend()
 fig.tight_layout();fig.savefig(ROOT/'results'/filename,dpi=170);plt.close(fig)
