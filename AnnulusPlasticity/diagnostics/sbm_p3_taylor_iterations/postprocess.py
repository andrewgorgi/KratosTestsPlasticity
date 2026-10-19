"""Compare untouched native runs and summarize replay strain diagnostics."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
summary=[]
for n in (20,40):
    baseline=json.loads((HERE/f'baseline/new/spans_{n}x{n}/result.json').read_text())
    audited=json.loads((HERE/f'audit/new/spans_{n}x{n}/result.json').read_text())
    assert baseline['status']==audited['status']
    assert baseline['metrics']==audited['metrics']
    folder=HERE/f'audit/audit_n{n}'
    steps=json.loads((folder/'steps.json').read_text())
    for path in sorted(folder.glob('iterations_t*.csv')):
        rows=list(csv.DictReader(path.open()));end=rows[-1];rr=rows[1:]
        R=np.array([float(r['residual_norm']) for r in rows])
        changes=np.diff(R);rising=changes>0
        ratios=np.array([float(r['median_second_term_over_first_order']) for r in rr])
        corr=float(np.corrcoef(np.log10(np.maximum(R[1:],1e-30)),ratios)[0,1]) if len(rr)>2 and np.std(ratios)>0 else None
        item=dict(mesh=n,time=float(end['time']),state=end['state'],iterations=int(end['iteration']),
            instrumented_status_and_metrics_identical=True,
            final_residual=float(R[-1]),peak_residual=float(R.max()),
            final_median_term2_ratio=float(end['median_second_term_over_first_order']),
            peak_median_term2_ratio=float(max(ratios)),
            residual_term2_correlation=corr,
            median_term2_on_residual_rises=float(np.median(ratios[rising])) if rising.any() else None,
            median_term2_on_residual_decreases=float(np.median(ratios[~rising])) if (~rising).any() else None,
            final_median_errors_q1_q2_q5=[float(end[f'median_error_q{k}_over_exact']) for k in (1,2,5)],
            final_median_q2_vs_local_continuation=float(end['median_q2_vs_local_continuation']),
            final_crossing_knots=int(end['crosses_knot']),final_crossing_exact_front=int(end['crosses_exact_front']),
            final_phase_mismatches=int(end['robust_phase_mismatch']),
            final_crossing_median_local_vs_global=float(end['crossing_median_local_vs_global']) if end['crossing_median_local_vs_global'] else None,
            final_same_span_median_local_vs_global=float(end['same_span_median_local_vs_global']) if end['same_span_median_local_vs_global'] else None,
            max_full_taylor_validation=float(max(float(r['max_full_taylor_validation']) for r in rows)))
        pairpath=path.with_name(path.name.replace('iterations_', 'pairs_').replace('.csv', '_iterations.csv'))
        pairs=list(csv.DictReader(pairpath.open()))
        last=[r for r in pairs if int(r['iteration'])==int(end['iteration'])]
        item['missing_global_support_pairs']=sum(int(r['missing_global_coefficients'])>0 for r in last)
        item['final_pairs_term2_larger_than_order1']=sum(float(r['second_term_over_first_order'])>1 for r in last)
        item['by_knot_group']={}
        for crossing in (False,True):
            group=[r for r in last if (r['crosses_knot']=='True')==crossing]
            item['by_knot_group']['crossing' if crossing else 'not_crossing']=dict(count=len(group),
                median_error_q1_q2_q5=[float(np.median([float(r[f'error_q{k}_over_exact']) for r in group])) for k in (1,2,5)],
                median_term2_ratio=float(np.median([float(r['second_term_over_first_order']) for r in group])))
        summary.append(item)
        if end['state']=='failed':
            fig,axes=plt.subplots(3,1,figsize=(9,8),sharex=True)
            x=[int(r['iteration']) for r in rows]
            axes[0].semilogy(x,R,label='Norma residuo globale');axes[0].legend()
            axes[1].plot(x,[float(r['median_second_term_over_first_order']) for r in rows],label='Mediana ||termine 2|| / ||Taylor 1||')
            axes[1].legend()
            for k in (1,2,5):axes[2].semilogy(x,[max(float(r[f'median_error_q{k}_over_exact']),1e-16) for r in rows],label=f'Errore strain mediano, Taylor {k}')
            axes[2].legend();axes[2].set_xlabel('Iterazione nativa dopo line search')
            for ax in axes:ax.grid(alpha=.25)
            fig.suptitle(f'SBM p3, {n}×{n}, passo fallito t={end["time"]}')
            fig.tight_layout();fig.savefig(HERE/f'iterations_n{n}.png',dpi=160);plt.close(fig)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
