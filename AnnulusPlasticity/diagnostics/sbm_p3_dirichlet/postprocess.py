"""Compare p3 Dirichlet cycles with the frozen mixed-boundary baseline."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'sbm_cyclic_old_vs_new_p3'
VARIANTS = {'free':(-1., .00025), 'free_dt125':(-1., .000125), 'penalty_100E':(6700000., .00025)}


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows, histories, checks = [], [], []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for ax, n in zip(axes, (20, 40)):
        baseline_case = BASE/'new'/f'spans_{n}x{n}'
        baseline = json.loads((baseline_case/'result.json').read_text())
        entries = [('mixed', baseline, -1., .00025)]
        for name, (penalty, dt) in VARIANTS.items():
            case = HERE/name/'new'/f'spans_{n}x{n}'
            data = json.loads((case/'result.json').read_text())
            entries.append((name, data, penalty, dt))
            actual = json.loads((case/'results'/'parameters_used.json').read_text())
            expected = json.loads((baseline_case/'results'/'parameters_used.json').read_text())
            replacements = 0
            for modeler in expected['modelers']:
                for layer in modeler['Parameters'].get('link_layer_to_condition_name', []):
                    if layer['condition_name']=='SbmLoadSolidCondition':
                        layer['condition_name']='SbmSolidCondition'
                        replacements += 1
            assert replacements==1
            expected['solver_settings']['time_stepping']['time_step']=dt
            for parameters in (actual, expected):
                parameters['solver_settings']['material_import_settings'].pop('materials_filename')
            mat = json.loads((case/'results'/'materials_used.json').read_text())
            expected_mat = json.loads((baseline_case/'results'/'materials_used.json').read_text())
            for prop in expected_mat['properties']:
                prop['Material']['Variables']['PENALTY_FACTOR']=penalty
            counts = json.loads((case/'results'/'conditions.json').read_text())
            no_loads = not any('SbmLoadSolidCondition' in key for key in counts)
            assert no_loads and any('SbmSolidCondition' in key for key in counts)
            with np.load(case/'results'/'solution.npz') as solution, np.load(baseline_case/'results'/'solution.npz') as old:
                same_points = np.array_equal(solution['points'], old['points'])
                same_weights = np.array_equal(solution['weights'], old['weights'])
            check = dict(variant=name, spans=n, parameters_match_expected=actual==expected,
                         materials_match_expected=mat==expected_mat, no_neumann_conditions=no_loads,
                         same_sampling_points=same_points, same_sampling_weights=same_weights,
                         conditions=counts)
            assert all(check[key] for key in ('parameters_match_expected', 'materials_match_expected',
                'no_neumann_conditions', 'same_sampling_points', 'same_sampling_weights')), check
            checks.append(check)
        for name, data, penalty, dt in entries:
            status, metrics = data['status'], data['metrics']
            row = dict(variant=name, spans=n, penalty=penalty, dt=dt, status=status['status'],
                       completed_steps=len(metrics), last_converged_time=metrics[-1]['time'],
                       failed_time=status['time'] if status['status']=='failed' else '',
                       L2_first_peak=next((m['L2_amplitude_relative'] for m in metrics if abs(m['time']-.001)<1e-10), ''),
                       L2_final=metrics[-1]['L2_amplitude_relative'] if status['status']=='converged' else '',
                       H1_final=metrics[-1]['H1_bulk_amplitude_relative'] if status['status']=='converged' else '',
                       max_L2=max(m['L2_amplitude_relative'] for m in metrics),
                       max_iterations=max(m['iterations'] for m in metrics))
            rows.append(row)
            histories += [dict(variant=name, spans=n, **m) for m in metrics]
            labels = {'mixed':'Mista, penalty-free', 'free':'Dirichlet, penalty-free',
                      'free_dt125':'Dirichlet, penalty-free, dt/2', 'penalty_100E':'Dirichlet, penalty 100E'}
            ax.semilogy([m['time']*1000 for m in metrics], [m['L2_amplitude_relative'] for m in metrics],
                        '.-', label=labels[name])
        ax.set(title=f'p3, mesh {n}×{n}', xlabel='Tempo [ms]', ylabel='L2 / norma al primo picco')
        ax.grid(True, which='both', alpha=.3)
        ax.legend(fontsize=8)
    fig.suptitle('Sostituzione della Neumann interna con Dirichlet analitica')
    fig.savefig(HERE/'comparison.png', dpi=180)
    plt.close(fig)
    write_csv(HERE/'summary.csv', rows)
    write_csv(HERE/'errors_by_time.csv', histories)
    metadata = json.loads((BASE/'build.json').read_text())
    source = Path(metadata['source'])
    library = source.parents[3]/'bin/Release/libs/libKratosIgaCore.so'
    diagnostic = metadata['libraries']['new']
    assert hashlib.sha256(Path(diagnostic['path']).read_bytes()).hexdigest()==diagnostic['sha256']
    hashes = dict(source_current=hashlib.sha256(source.read_bytes()).hexdigest(),
                  source_frozen=metadata['sha256'],
                  installed_library_current=hashlib.sha256(library.read_bytes()).hexdigest(),
                  installed_library_frozen=metadata['installed_library_sha256'],
                  diagnostic_library=diagnostic)
    (HERE/'checks.json').write_text(json.dumps(dict(cases=checks, hashes=hashes,
        source_matches_frozen=hashes['source_current']==hashes['source_frozen'],
        installed_library_matches_frozen=hashes['installed_library_current']==hashes['installed_library_frozen'],
        diagnostic_library_matches_frozen=True), indent=2)+'\n')
    print(json.dumps(rows, indent=2))


if __name__=='__main__':
    main()
