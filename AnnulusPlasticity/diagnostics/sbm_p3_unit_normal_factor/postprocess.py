"""Summarize and plot accepted steps of original versus unit-factor SBM runs."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
summary=[]
fig,axes=plt.subplots(2,2,figsize=(11,7),sharex='col')
for col,n in enumerate((20,40)):
    for mode,label in [('old','Prodotto delle normali'),('new','Fattore 1')]:
        result=json.loads((HERE/mode/f'spans_{n}x{n}/result.json').read_text());status=result['status'];metrics=result['metrics']
        row=dict(mesh=n,variant=label,status=status['status'],stop_time=status.get('time'),completed_steps=status.get('completed_steps'),
            last_accepted_time=metrics[-1]['time'] if metrics else None,
            last_L2=metrics[-1]['L2_amplitude_relative'] if metrics else None,
            max_L2=max((r['L2_amplitude_relative'] for r in metrics),default=None),
            max_NR_accepted=max((r['iterations'] for r in metrics),default=None))
        summary.append(row)
        t=[r['time'] for r in metrics]
        line,=axes[0,col].semilogy(t,[r['L2_amplitude_relative'] for r in metrics],'.-',label=label)
        axes[1,col].plot(t,[r['iterations'] for r in metrics],'.-',color=line.get_color(),label=label)
        if status['status']=='failed':axes[0,col].axvline(status['time'],ls=':',color=line.get_color(),alpha=.7)
    axes[0,col].set_title(f'SBM p3, {n}×{n}')
    axes[0,col].set_ylabel('Errore L2 / ampiezza analitica');axes[1,col].set_ylabel('Iterazioni NR dei passi accettati');axes[1,col].set_xlabel('Tempo')
    for ax in axes[:,col]:ax.grid(alpha=.25);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(HERE/'comparison.png',dpi=160)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (HERE/'summary.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=summary[0]);writer.writeheader();writer.writerows(summary)
# Compare errors only at common accepted times; a failed state's error is never evaluated.
common=[]
for n in (20,40):
    data={mode:json.loads((HERE/mode/f'spans_{n}x{n}/result.json').read_text()) for mode in ('old','new')}
    old={round(r['time'],12):r for r in data['old']['metrics']}
    for r in data['new']['metrics']:
        if round(r['time'],12) in old:
            prev=old[round(r['time'],12)]
            common.append(dict(mesh=n,time=r['time'],L2_dot=prev['L2_amplitude_relative'],L2_unit=r['L2_amplitude_relative'],iterations_dot=prev['iterations'],iterations_unit=r['iterations']))
(HERE/'common_steps.json').write_text(json.dumps(common,indent=2)+'\n')
print(json.dumps(summary,indent=2));print('COMMON_STEPS',json.dumps(common))
