#!/usr/bin/env python3
"""Plot numerical historical-window sensitivity. No constitutive hysteresis inferred."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'results'
with (R/'event_echo_windows.csv').open() as f: rows=list(csv.DictReader(f))
rows=[r for r in rows if r['tolerance']=='standard']
hs=[.25,.5,1.,2.]
fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
for event in range(1,5):
    d=sorted((r for r in rows if int(r['event'])==event),key=lambda r:float(r['window_xi']))
    ax.plot(hs,[float(r['direction_line_angle_deg']) for r in d],marker='o',label=f'Event {event} — xi {float(d[0]["xi"]):.3f}')
ax.set(xlabel='Retained recent history h (xi)',ylabel='Line angle from complete history (degrees)',title='Q550 v18.14 — Effect of removing older source history')
ax.set_xticks(hs);ax.set_ylim(bottom=0);ax.grid(alpha=.25);ax.legend()
fig.savefig(R/'echo_angle_by_window.png',dpi=175);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
for event in range(1,5):
    d=sorted((r for r in rows if int(r['event'])==event),key=lambda r:float(r['window_xi']))
    ax.plot(hs,[float(r['old_exclusive_over_full_radius']) for r in d],marker='o',label=f'Event {event}')
ax.set(xlabel='Retained recent history h (xi)',ylabel='Old direction outside recent plane / full radius',title='Q550 v18.14 — Exclusive old-history directional contribution')
ax.set_xticks(hs);ax.grid(alpha=.25);ax.legend()
fig.savefig(R/'old_history_exclusive_echo.png',dpi=175);plt.close(fig)
with (R/'history_window_profiles.csv').open() as f: prof=list(csv.DictReader(f))
prof=[r for r in prof if r['tolerance']=='standard']
fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
for h in hs:
    rr=[r for r in prof if abs(float(r['window_xi'])-h)<1e-8]
    ax.plot([float(r['xi']) for r in rr],[float(r['line_angle_deg']) for r in rr],label=f'last {h:g} xi',lw=1.3)
with (R/'summary.json').open() as f: info=json.load(f)
for q in info['by_tolerance'][0]['events_xi']:
    ax.axvline(q,ls=':',alpha=.3)
ax.set(xlabel='xi',ylabel='Line angle from complete-history plane (deg)',title='Q550 v18.14 — Source history sensitivity over the trajectory')
ax.grid(alpha=.2);ax.legend()
fig.savefig(R/'history_echo_profile.png',dpi=175);plt.close(fig)
print('Plots saved:',*[p.name for p in R.glob('*.png')])
