"""Simple diagnostic figures from real Q550 frozen checkpoint results."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def read(n):
    with (ROOT/'results'/n).open() as f:return list(csv.DictReader(f))

w=read('direct_time_window_observability.csv')
fig,ax=plt.subplots(figsize=(8.4,4.8))
for obs in ['material3','raw4','raw6']:
    rr=sorted([r for r in w if r['observer']==obs],key=lambda x:float(x['half_window']))
    ax.semilogy([float(r['half_window']) for r in rr],[float(r['condition']) for r in rr],marker='o',label=obs)
ax.set_xlabel('Observation half-window (xi)')
ax.set_ylabel('Condition number (nine direct time samples)')
ax.set_title('Q550 v18.21 | Hidden-state visibility versus observation length')
ax.legend();ax.grid(alpha=.2);fig.tight_layout()
fig.savefig(ROOT/'results/direct_window_conditioning.png',dpi=175);plt.close(fig)

c=read('jet_noise_singular_spectrum.csv')
fig,ax=plt.subplots(figsize=(8.4,4.8))
for obs in ['material3','raw4','raw6']:
    rr=sorted([r for r in c if r['observer']==obs and float(r['xi'])==4.],key=lambda x:float(x['tau']))
    ax.loglog([float(r['tau']) for r in rr],[float(r['condition']) for r in rr],marker='o',label=obs)
ax.set_xlabel('Jet time scaling tau (xi)')
ax.set_ylabel('Naive scaled-jet condition number')
ax.set_title('Q550 v18.21 | Changing derivative units changes naive conditioning')
ax.legend();ax.grid(alpha=.2);fig.tight_layout()
fig.savefig(ROOT/'results/jet_scale_naive_conditioning.png',dpi=175);plt.close(fig)

s=read('noise_recovery_trials.csv')
fig,ax=plt.subplots(figsize=(8.4,4.8))
for obs in ['material3','raw4','raw6']:
    rr=sorted([r for r in s if r['observer']==obs and float(r['tau'])==.1],key=lambda x:float(x['noise_sigma']))
    ax.loglog([float(r['noise_sigma']) for r in rr],[float(r['empirical_rmse']) for r in rr],marker='o',label=obs)
ax.set_xlabel('Assumed independent noise std per scaled jet entry')
ax.set_ylabel('RMS error reconstructing weakest unit 10D state direction')
ax.set_title('Q550 v18.21 | Synthetic noise, not physical instrument data')
ax.legend();ax.grid(alpha=.2);fig.tight_layout()
fig.savefig(ROOT/'results/synthetic_noise_recovery.png',dpi=175);plt.close(fig)
print('three diagnostic PNG files generated')
