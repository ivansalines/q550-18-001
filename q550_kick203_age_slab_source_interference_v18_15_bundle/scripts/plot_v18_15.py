#!/usr/bin/env python3
"""Render only computed gauge-declared v18.15 results."""
import csv,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'results'
def read(name):
    with (P/name).open() as f:return list(csv.DictReader(f))
def main():
    x=read('event_age_slabs.csv')
    # Four separate plots, never merge distinct diagnostics in a subplot.
    for event in range(1,5):
        r=[a for a in x if a['tolerance']=='standard' and int(a['event'])==event and abs(float(a['width'])-.25)<1e-9]
        age=np.array([float(v['age_mid']) for v in r]); val=np.array([float(v['signed_axial_over_full']) for v in r]);
        fig,ax=plt.subplots(figsize=(9,3.7));ax.bar(age,val,width=[float(v['age_end'])-float(v['age_start']) for v in r])
        ax.axhline(0,linewidth=.8);ax.axvline(2,linestyle='--',linewidth=1)
        ax.set(xlabel='Source age (xi): recent at left, remote at right',ylabel='Signed slab contribution / final radius',title=f'Q550 v18.15 — event {event}: source-age echo, slabs 0.25 xi')
        fig.tight_layout();fig.savefig(P/f'age_slab_event_{event}.png',dpi=160);plt.close(fig)
    rr=read('resolution_convergence.csv')
    fig,ax=plt.subplots(figsize=(8,4))
    for i in range(1,5):
        a=sorted([z for z in rr if z['tolerance']=='standard' and int(z['event'])==i],key=lambda z:-float(z['width']))
        ax.plot([float(z['width']) for z in a],[float(z['destructive_fraction']) for z in a],marker='o',label=f'Event {i}')
    ax.set(xlabel='Source-age slab width (xi)',ylabel='Signed axial cancellation fraction',title='Resolution stress: how cancellation changes when history is subdivided')
    ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(P/'cancellation_resolution.png',dpi=160);plt.close(fig)
if __name__=='__main__':main()
