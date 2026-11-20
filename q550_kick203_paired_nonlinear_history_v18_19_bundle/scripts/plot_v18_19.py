"""Plot server nonlinear trajectories; no plot if server audit not yet run."""
from pathlib import Path
import csv,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results/paired_trajectories.csv'
with p.open() as f:rows=list(csv.DictReader(f))
fig,ax=plt.subplots(figsize=(9.5,5.5))
for name in ('material3_jet_a','material3_state_only_a','raw4_jet_a','raw4_state_only_a','raw6_state_only_a'):
 r=[x for x in rows if x['label']==name and float(x['epsilon'])==3e-6 and float(x['xi'])>=4.]
 if r:
  ax.plot([float(x['xi']) for x in r],[float(x['scaled_reduced_separation']) for x in r],label=name)
ax.set(xlabel='xi',ylabel='Reduced future separation / (2 epsilon)',
       title='Q550 v18.19 — Paired nonlinear futures at matched reduced present')
ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout()
fig.savefig(ROOT/'results/paired_future_discrimination.png',dpi=160);plt.close(fig)
