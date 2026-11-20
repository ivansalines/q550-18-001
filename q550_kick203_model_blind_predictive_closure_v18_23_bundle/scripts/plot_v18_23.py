"""Plot authentic frozen / optional original nonlinear experiment summaries."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
for scope in ('frozen','original'):
    file=ROOT/'results'/f'{scope}_model_comparison.csv'
    if not file.exists():continue
    rows=list(csv.DictReader(file.open()))
    plt.figure(figsize=(9,5))
    for observer in ('material3','raw4','raw6'):
        part=[r for r in rows if r['observer']==observer and float(r['noise_sigma'])==1e-4]
        vals={r['model']:float(r['test_rmse']) for r in part}
        if not vals:continue
        plt.plot(range(4),[vals[m] for m in ('instant','recent_0.10','history_0.50','per_trial_shuffled')],marker='o',label=observer)
    plt.xticks(range(4),['Instant','Last 0.10 ξ','Last 0.50 ξ','Shuffled history'])
    plt.ylabel('Held-out future RMSE (normalized diagnostic output)')
    plt.yscale('log');plt.grid(alpha=.2);plt.legend();plt.title(f'Q550 v18.23: {scope} sensor-only prediction')
    plt.tight_layout()
    plt.savefig(ROOT/'results'/f'{scope}_blind_prediction.png',dpi=150)
    plt.close()
