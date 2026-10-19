"""Summarize native branch switches and the actual Kratos utility experiments."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
stats=[];counts=[]
for n in (20,40):
    rows=list(csv.DictReader((HERE/f'local/points_n{n}.csv').open()))
    for iteration in (0,2):
        for group in ('bulk','surrogate','true'):
            rr=[r for r in rows if int(r['iteration'])==iteration and r['group']==group and float(r['eta'])==1e-8 and float(r['direction_response_norm'])>1e-8]
            counts.append(dict(mesh=n,iteration=iteration,group=group,points=len(rr),near_yield=sum(r['near_yield']=='True' for r in rr),
                switches_forward=sum(r['switches_forward']=='True' for r in rr),switches_backward=sum(r['switches_backward']=='True' for r in rr)))
            categories={'all':rr,'near_yield':[r for r in rr if r['near_yield']=='True'],
                'away_same_branch':[r for r in rr if r['near_yield']=='False' and r['switches_forward']=='False'],
                'changed_forward':[r for r in rr if r['switches_forward']=='True']}
            for label,selected in categories.items():
                if not selected:continue
                row=dict(mesh=n,iteration=iteration,group=group,category=label,points=len(selected))
                for key in ('analytic_forward_error','analytic_backward_error','utility_forward_order1_error','utility_forward_order2_error','utility_central_small_h_error'):
                    vals=np.array([float(r[key]) for r in selected if r[key]!=''])
                    weights=np.array([float(r['direction_response_norm']) for r in selected if r[key]!=''])
                    row[key+'_median']=float(np.median(vals)) if len(vals) else None
                    row[key+'_max']=float(vals.max()) if len(vals) else None
                    row[key+'_weighted_rms']=float(np.linalg.norm(vals*weights)/np.linalg.norm(weights)) if len(vals) else None
                stats.append(row)
(HERE/'local_counts.json').write_text(json.dumps(counts,indent=2)+'\n')
(HERE/'local_statistics.json').write_text(json.dumps(stats,indent=2)+'\n')
with (HERE/'local_counts.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=counts[0]);w.writeheader();w.writerows(counts)
results=[]
cases=[('utility_central',20),('utility_central',40),('utility_forward',20),('utility_central_h1e10',20),('utility_central_full',20)]
fig,axes=plt.subplots(1,2,figsize=(11,4.5))
for variant,n in cases:
    p=HERE/variant/'new'/f'spans_{n}x{n}/result.json'
    if not p.exists():raise RuntimeError(f'Pending result: {p}')
    d=json.loads(p.read_text());metrics=d['metrics']
    results.append(dict(variant=variant,mesh=n,**d['status'],last_L2=metrics[-1]['L2_amplitude_relative'] if metrics else None,
        max_L2=max((r['L2_amplitude_relative'] for r in metrics),default=None),max_NR=max((r['iterations'] for r in metrics),default=None),metrics=metrics))
    if variant in ('utility_central','utility_central_full'):
        label=f'{variant}, {n}×{n}'
        axes[0].semilogy([r['time'] for r in metrics],[r['L2_amplitude_relative'] for r in metrics],'.-',label=label)
        axes[1].plot([r['time'] for r in metrics],[r['iterations'] for r in metrics],'.-',label=label)
for ax in axes:ax.set_xlabel('Tempo');ax.grid(alpha=.25);ax.legend(fontsize=7)
axes[0].set_ylabel('L2 / ampiezza analitica');axes[1].set_ylabel('Iterazioni NR dei passi accettati')
fig.tight_layout();fig.savefig(HERE/'utility_runs.png',dpi=160)
(HERE/'solver_summary.json').write_text(json.dumps(results,indent=2)+'\n')
for r in results:print(json.dumps({k:v for k,v in r.items() if k!='metrics'}))
print('COUNTS',json.dumps(counts))
