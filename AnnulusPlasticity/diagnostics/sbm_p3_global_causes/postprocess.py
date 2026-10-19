"""Collect global-cause diagnostics and visualize the negative boundary mode."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'sbm_neumann_current'


def main():
    rows = []
    for path in sorted(HERE.glob('*/new/*/result.json')):
        data = json.loads(path.read_text())
        metrics, status = data['metrics'], data['status']
        rows.append(dict(variant=path.parts[-4], spans=data['spans'], status=status['status'],
                         failed_time=status['time'] if status['status']=='failed' else '',
                         completed_steps=len(metrics), last_converged_time=metrics[-1]['time'],
                         last_L2=metrics[-1]['L2_amplitude_relative'],
                         maximum_L2=max(m['L2_amplitude_relative'] for m in metrics)))
    with (HERE/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    audit = HERE/'diagnostic'/'audit_n20'
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    angle = np.linspace(0., 2*np.pi, 500)
    for ax, name in zip(axes, ('elastic', 'plastic')):
        summary = json.loads((audit/name/'summary.json').read_text())
        with np.load(audit/name/'mode.npz') as data:
            points, indices = np.unique(data['coordinates'], axis=0, return_inverse=True)
            amplitude = np.zeros(len(points)); np.add.at(amplitude, indices, data['mode']**2)
            amplitude = np.sqrt(amplitude)
        dots = ax.scatter(points[:, 0], points[:, 1], c=amplitude/amplitude.max(), s=25, cmap='magma', vmin=0., vmax=1.)
        for radius in (1., 2.):
            ax.plot(radius*np.cos(angle), radius*np.sin(angle), color='black', lw=.7)
        ax.set(aspect='equal', title=f"t={summary['time']:g}: λ_sym,min={summary['symmetric_min_eigenvalue']:.0f}", xlabel='x', ylabel='y')
        figure.colorbar(dots, ax=ax, label='Ampiezza del modo / massimo')
    figure.suptitle('p3, 20×20: modo più negativo della parte simmetrica della tangente')
    figure.savefig(HERE/'negative_mode.png', dpi=180); plt.close(figure)
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    restart = json.loads((audit/'residual_armijo.json').read_text())
    values = [r['residual_norm'] for r in restart['iterations']]+[restart['final_residual_norm']]
    axes[0].semilogy(range(len(values)), values, '.-', label='20×20, ripartenza a t=0.00075')
    stalled = json.loads((HERE/'residual_ls'/'newton_trace_n40'/'t0.00075000.json').read_text())
    axes[0].semilogy([r['iteration']-1 for r in stalled['iterations']], [r['residual_norm'] for r in stalled['iterations']], '.-', label='40×40, passo t=0.00075')
    axes[0].set(xlabel='Aggiornamento Newton', ylabel='Norma del residuo', title='Armijo sulla norma del residuo')
    axes[0].legend(fontsize=8)
    directional = json.loads((HERE/'residual_ls'/'newton_trace_n20'/'directional_t0.00100000.json').read_text())
    axes[1].semilogx([r['alpha'] for r in directional], [r['relative_directional_derivative_mismatch'] for r in directional], 'o-')
    axes[1].set(xlabel='Frazione della correzione Newton', ylabel='Scarto relativo della derivata direzionale',
                title='20×20, inizio passo t=0.001')
    for ax in axes: ax.grid(True, which='both', alpha=.3)
    figure.savefig(HERE/'nonlinear_controls.png', dpi=180); plt.close(figure)
    metadata = json.loads((BASE/'build.json').read_text())
    source = Path(metadata['source'])
    library = Path(metadata['libraries']['new']['path'])
    checks = dict(source_matches_frozen=hashlib.sha256(source.read_bytes()).hexdigest()==metadata['sha256'],
                  installed_library_matches_frozen=hashlib.sha256(library.read_bytes()).hexdigest()==metadata['libraries']['new']['sha256'])
    baseline = json.loads((BASE/'new'/'spans_20x20'/'result.json').read_text())
    diagnostic = json.loads((HERE/'diagnostic'/'new'/'spans_20x20'/'result.json').read_text())
    checks['global_audit_preserves_native_status'] = baseline['status']==diagnostic['status']
    checks['global_audit_preserves_native_metrics'] = baseline['metrics']==diagnostic['metrics']
    derivatives = json.loads((HERE/'derivatives10'/'new'/'spans_20x20'/'result.json').read_text())
    checks['derivatives10_preserves_native_status'] = baseline['status']==derivatives['status']
    checks['derivatives10_preserves_native_metrics'] = baseline['metrics']==derivatives['metrics']
    for name in ('diagnostic', 'newton', 'limited_ls', 'lambda05', 'lambda1', 'quadrature8', 'residual_ls', 'derivatives10'):
        assert (HERE/name/'materials_input.json').read_bytes()==(BASE/'materials_input.json').read_bytes()
    checks['all_materials_unchanged'] = True
    assert all(checks.values()), checks
    (HERE/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps(dict(configurations=rows, checks=checks), indent=2))


if __name__=='__main__': main()
