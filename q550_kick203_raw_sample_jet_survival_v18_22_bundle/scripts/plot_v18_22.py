"""Draw plots from actual frozen Q550 tangent Monte Carlo CSVs."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'results/raw_sample_reconstruction_trials.csv').open() as f:rows=list(csv.DictReader(f))
fig,axes=plt.subplots(1,3,figsize=(13.5,3.6),sharey=True)
for ax,label in zip(axes,('material3_jet2','raw4_jet1','raw6_jet1')):
    for deg in [3,5]:
        r=[v for v in rows if v['label']==label and float(v['noise_sigma'])==.0001 and int(v['degree'])==deg]
        r.sort(key=lambda z:float(z['window_halfwidth']))
        x=[float(z['window_halfwidth']) for z in r]
        if deg==3:ax.plot(x,[float(z['direct_rms']) for z in r],marker='o',label='Raw window (direct)')
        ax.plot(x,[float(z['jet_rms']) for z in r],marker='o',label=f'Polynomial jet (degree {deg})')
    ax.set_title(label)
    ax.set_xlabel('Half-window in xi')
    ax.grid(alpha=.23)
    ax.set_yscale('log')
axes[0].set_ylabel('Full 10D reconstruction RMS (unit seed)')
axes[-1].legend(loc='best',fontsize=8)
fig.suptitle('Q550 v18.22 — equal raw-sample noise, sigma=1e-4; 33 samples')
fig.tight_layout();fig.savefig(ROOT/'results/raw_vs_jet_noise_survival.png',dpi=155);plt.close(fig)
