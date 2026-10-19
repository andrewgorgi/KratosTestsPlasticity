"""Summarize one-sided derivatives and fixed-history projected restarts."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
summary=[]
for n in (20,40):
    folder=HERE/f'audit_n{n}'
    baseline=json.loads((folder/'baseline.json').read_text());warm=json.loads((folder/'warm_start.json').read_text())
    directions=json.loads((folder/'native_directions.json').read_text())
    rows=list(csv.DictReader((folder/'directional.csv').open()))
    indices=sorted(set(int(r['iteration']) for r in rows))
    details=[]
    fig,axes=plt.subplots(2,len(indices),figsize=(13,7),sharex=True)
    for col,index in enumerate(indices):
        for group in ('bulk','dirichlet','neumann','total'):
            selected=[r for r in rows if int(r['iteration'])==index and r['group']==group]
            line,=axes[0,col].loglog([float(r['eta']) for r in selected],[float(r['forward_relative_mismatch']) for r in selected],'.-',label=group)
            axes[1,col].loglog([float(r['eta']) for r in selected],[float(r['backward_relative_mismatch']) for r in selected],'.-',color=line.get_color(),label=group)
            best=next(r for r in selected if float(r['eta'])==1e-6)
            details.append(dict(iteration=index,group=group,eta=1e-6,forward_mismatch=float(best['forward_relative_mismatch']),backward_mismatch=float(best['backward_relative_mismatch']),
                predicted_phi_slope=float(best['predicted_phi_slope']),forward_phi_slope=float(best['forward_phi_slope'])))
        axes[0,col].set_title(f'Iterato {index}');axes[0,col].set_ylabel('Errore derivata in avanti');axes[1,col].set_ylabel('Errore derivata indietro');axes[1,col].set_xlabel('eta × aggiornamento effettivo Newton')
        for ax in axes[:,col]:ax.grid(alpha=.25);ax.legend(fontsize=7)
    fig.suptitle(f'SBM p3 {n}×{n}, t=0.00075 — storia precedente invariata')
    fig.tight_layout();fig.savefig(HERE/f'directional_n{n}.png',dpi=160);plt.close(fig)
    result=json.loads((HERE/f'new/spans_{n}x{n}/result.json').read_text())
    previous=json.loads((HERE.parent/f'sbm_p3_taylor_iterations/baseline/new/spans_{n}x{n}/result.json').read_text())
    accepted_before=[r for r in result['metrics'] if r['time']<.00075-1e-12]
    assert accepted_before==previous['metrics']
    initial=[r for r in rows if r['iteration']=='0' and float(r['eta'])==1e-6]
    norms={r['group']:float(r['predicted_derivative_norm']) for r in initial}
    cancellation=sum(norms[k] for k in ('bulk','dirichlet','neumann'))/norms['total']
    item=dict(initial_direction_cancellation_ratio=cancellation,mesh=n,baseline=baseline,warm_start=warm,native_directions=directions,derivatives=details,previous_accepted_steps_identical=True)
    fresh=json.loads((HERE/f'fresh/audit_n{n}/warm_start.json').read_text())
    fresh_result=json.loads((HERE/f'fresh/new/spans_{n}x{n}/result.json').read_text())
    assert [r for r in fresh_result['metrics'] if r['time']<.00075-1e-12]==previous['metrics']
    item['fresh_warm_start']=fresh
    item['fresh_restart_matches_exactly']=warm==fresh
    summary.append(item)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axes=plt.subplots(1,2,figsize=(10,4))
for n in (20,40):
    scan=list(csv.DictReader((HERE/f'audit_n{n}/projection_scan.csv').open()))
    x=[float(r['fraction']) for r in scan]
    axes[0].semilogy(x,[float(r['residual_norm']) for r in scan],'.-',label=f'{n}×{n}')
    axes[1].semilogy(x,[float(r['L2_relative_current']) for r in scan],'.-',label=f'{n}×{n}')
for ax in axes:ax.set_xlabel('Frazione verso la proiezione L2 analitica');ax.legend();ax.grid(alpha=.25)
axes[0].set_ylabel('Norma residuo');axes[1].set_ylabel('Errore L2 relativo al campo corrente')
fig.tight_layout();fig.savefig(HERE/'projection_scan.png',dpi=160)
for item in summary:
    print('MESH',item['mesh']);print('BASELINE',json.dumps(item['baseline']));print('WARM',json.dumps(item['warm_start']));print('NATIVE',json.dumps(item['native_directions']));print('DERIVATIVES',json.dumps(item['derivatives']))
