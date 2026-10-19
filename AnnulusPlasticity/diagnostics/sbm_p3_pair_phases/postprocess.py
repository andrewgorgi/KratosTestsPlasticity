"""Summarize robust phase mismatches in native Newton iteration snapshots."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
summary = []
for n, ax in zip((20, 40), axes):
    rows = list(csv.DictReader((HERE/f'audit_n{n}/iterations_t0.00075000_failed.csv').open()))[1:]
    x = [int(r['iteration']) for r in rows]
    a = [int(r['robust_elastic_s_plastic_t']) for r in rows]
    b = [int(r['robust_plastic_s_elastic_t']) for r in rows]
    ax.plot(x, a, '.-', label='Surrogate elastico / true plastico', linewidth=1)
    ax.plot(x, b, '.-', label='Surrogate plastico / true elastico', linewidth=1)
    ax.set_title(f'{n}×{n}, p3 — passo fallito t = 0.00075')
    ax.set_ylabel('Coppie discordanti')
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    current=json.loads((HERE/f'new/spans_{n}x{n}/result.json').read_text())
    baseline=json.loads((HERE.parent/f'sbm_neumann_current/new/spans_{n}x{n}/result.json').read_text())
    assert current['status']==baseline['status']
    assert current['metrics']==baseline['metrics']
    summary.append(dict(spans=n, points=int(rows[0]['points']),
        failed_time=float(rows[0]['time']), iterations_with_mismatch=sum(v+w>0 for v,w in zip(a,b)),
        max_elastic_s_plastic_t=max(a), max_plastic_s_elastic_t=max(b),
        final_elastic_s_plastic_t=a[-1], final_plastic_s_elastic_t=b[-1],
        native_status_and_metrics_unchanged=True))
axes[-1].set_xlabel('Iterazione Newton (dopo aggiornamento e line search)')
fig.tight_layout()
fig.savefig(HERE/'phase_pairs.png', dpi=180)
(HERE/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
