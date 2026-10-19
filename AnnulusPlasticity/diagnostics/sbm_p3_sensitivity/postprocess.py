"""Summarize time-step and penalty experiments, including accepted-state errors."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigvals, eigvalsh
from scipy.sparse import load_npz

HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'sbm_cyclic_old_vs_new_p3'


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    cases, checks = [], []
    paths = [(BASE/'new'/f'spans_{n}x{n}'/'result.json', 'baseline', .00025, -1.) for n in (20, 40)]
    for folder in sorted(HERE.iterdir()):
        if not (folder/'materials_input.json').exists() or folder.name == 'jacobian_check':
            continue
        parameters = json.loads((folder/'parameters_input.json').read_text())
        material = json.loads((folder/'materials_input.json').read_text())
        dt = parameters['solver_settings']['time_stepping']['time_step']
        penalty = material['properties'][0]['Material']['Variables']['PENALTY_FACTOR']
        paths += [(p, folder.name, dt, penalty) for p in folder.glob('new/spans_*/result.json')]
    for path, name, dt, penalty in paths:
        data = json.loads(path.read_text())
        metrics, status = data['metrics'], data['status']
        row = dict(variant=name, spans=data['spans'], dt=dt, penalty=penalty,
                   status=status['status'], completed_steps=len(metrics),
                   last_converged_time=metrics[-1]['time'] if metrics else 0.,
                   failed_time=status['time'] if status['status']=='failed' else '',
                   L2_at_first_peak=next((m['L2_amplitude_relative'] for m in metrics if abs(m['time']-.001)<1e-10), ''),
                   last_L2=metrics[-1]['L2_amplitude_relative'] if metrics else '',
                   max_L2=max((m['L2_amplitude_relative'] for m in metrics), default=0.),
                   max_iterations=max((m['iterations'] for m in metrics), default=0))
        cases.append((row, metrics))
        actual = json.loads((path.parent/'results'/'parameters_used.json').read_text())
        baseline = json.loads((BASE/'new'/f"spans_{data['spans']}x{data['spans']}"/'results'/'parameters_used.json').read_text())
        for parameters in (actual, baseline):
            parameters['solver_settings']['time_stepping'].pop('time_step')
            parameters['solver_settings']['material_import_settings'].pop('materials_filename')
        mat = json.loads((path.parent/'results'/'materials_used.json').read_text())
        baseline_mat = json.loads((BASE/'new'/f"spans_{data['spans']}x{data['spans']}"/'results'/'materials_used.json').read_text())
        for material in (mat, baseline_mat):
            for prop in material['properties']:
                prop['Material']['Variables'].pop('PENALTY_FACTOR')
        check = dict(variant=name, spans=data['spans'], parameters_only_dt_changed=actual==baseline,
                     material_only_penalty_changed=mat==baseline_mat)
        assert check['parameters_only_dt_changed'] and check['material_only_penalty_changed'], check
        checks.append(check)
    rows = [r for r, _ in cases]
    write_csv(HERE/'summary.csv', rows)
    histories = [dict(variant=row['variant'], spans=row['spans'], dt=row['dt'], penalty=row['penalty'], **m)
                 for row, metrics in cases for m in metrics]
    write_csv(HERE/'errors_by_time.csv', histories)
    metadata = json.loads((BASE/'build.json').read_text())
    source = Path(metadata['source'])
    library = source.parents[3]/'bin/Release/libs/libKratosIgaCore.so'
    assert hashlib.sha256(source.read_bytes()).hexdigest()==metadata['sha256']
    assert hashlib.sha256(library.read_bytes()).hexdigest()==metadata['installed_library_sha256']
    (HERE/'checks.json').write_text(json.dumps(dict(cases=checks, source_unchanged=True,
                                                  installed_library_unchanged=True), indent=2)+'\n')
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    predicates = [lambda r:r['spans']==20 and r['penalty']==-1,
                  lambda r:r['spans']==20 and r['dt']==.00025,
                  lambda r:r['spans']==20 and r['penalty']>0 and r['dt']<.00025,
                  lambda r:r['spans']==40]
    titles = ['20×20: time step, penalty-free', '20×20: penalty, dt=0.00025',
              '20×20: penalty e passo ridotto', '40×40: controllo su mesh più fine']
    for ax, predicate, title in zip(axes.flat, predicates, titles):
        for row, metrics in cases:
            if predicate(row):
                visible = [m for m in metrics if m['time'] <= .001+1e-12]
                ax.semilogy([m['time']*1000 for m in visible],
                            [m['L2_amplitude_relative'] for m in visible], '.-', label=row['variant'])
        ax.set(title=title, xlabel='Tempo [ms]', ylabel='L2 / norma al primo picco', xlim=(0, 1.02))
        ax.grid(True, which='both', alpha=.3)
        ax.legend(fontsize=7)
    fig.suptitle('SBM p3, legge al punto proiettato: errori dei soli passi accettati')
    fig.savefig(HERE/'sensitivity.png', dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 9), constrained_layout=True)
    for i, (row, metrics) in enumerate(cases):
        ax.barh(i, row['last_converged_time']*1000, color='tab:blue', alpha=.65)
        if row['failed_time'] != '':
            ax.plot(row['failed_time']*1000, i, 'rx')
    ax.set_yticks(range(len(rows)), [f"{r['variant']} — {r['spans']}×{r['spans']}" for r in rows])
    ax.set(xlabel='Tempo [ms]', xlim=(0, 6.2), title='Ultimo passo accettato; × rosso = arresto Newton')
    ax.axvline(6., color='black', linestyle=':', label='Fine ciclo previsto')
    ax.grid(True, axis='x', alpha=.3)
    ax.legend()
    fig.savefig(HERE/'availability.png', dpi=180)
    plt.close(fig)
    root = HERE/'jacobian_check'
    spectra = {}
    for state in ('elastic_state', 'converged_state', 'failed_state'):
        matrix = load_npz(root/state/'cycle_matrix.npz').toarray()
        ev = eigvals(matrix)
        symmetric = eigvalsh((matrix+matrix.T)/2)
        threshold = 1e-9*np.max(np.abs(ev))
        spectra[state] = dict(min_real_eigenvalue=float(ev.real.min()),
            max_real_eigenvalue=float(ev.real.max()), negative_real_eigenvalues=int((ev.real < -threshold).sum()),
            negative_threshold=float(threshold), symmetric_part_min_eigenvalue=float(symmetric.min()),
            non_symmetry_relative=float(np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix)))
        if state == 'elastic_state':
            zero = load_npz(root/state/'elastic_matrix.npz').toarray()
            spectra[state]['relative_difference_from_zero_displacement_matrix'] = float(
                np.linalg.norm(matrix-zero)/np.linalg.norm(zero))
    (root/'spectra.json').write_text(json.dumps(spectra, indent=2)+'\n')
    print(json.dumps(dict(cases=len(cases), completed_cycles=sum(r['status']=='converged' for r in rows),
                         spectra=spectra), indent=2))


if __name__ == '__main__':
    main()
