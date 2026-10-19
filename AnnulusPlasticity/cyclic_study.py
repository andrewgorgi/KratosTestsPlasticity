"""Cyclic annular torsion: convergence and contour histories for SBM/Gap-SBM."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import annulus_study as base
import cyclic_analytic as exact

ROOT, KM = base.ROOT, base.KM


def errors(samples):
    """Fixed first-peak normalization, including at zero torque/residual states."""
    points, w = samples['points'], samples['weights']
    u, grad, _, _ = exact.evaluate(points, samples['time'])
    u0, grad0, _, _ = exact.evaluate(points, exact.TIMES[1])
    bulk = np.isfinite(samples['gradient']).all(axis=1)
    l2 = np.sum(w[:, None]*(samples['displacement']-u)**2)
    l2norm = np.sum(w[:, None]*u0**2)
    h1 = np.sum(w[bulk, None]*((samples['displacement'][bulk]-u[bulk])**2))
    h1 += np.sum(w[bulk, None]*(samples['gradient'][bulk]-grad[bulk])**2)
    h1norm = np.sum(w[bulk, None]*u0[bulk]**2)+np.sum(w[bulk, None]*grad0[bulk]**2)
    if not np.isfinite([l2, h1, l2norm, h1norm]).all() or min(l2, h1) < 0:
        raise RuntimeError('Invalid squared cyclic error (nonfinite values or signed gap quadrature).')
    if min(l2norm, h1norm) <= 0:
        raise RuntimeError('Zero first-peak normalization.')
    return dict(time=float(samples['time']),
        load_factor=float(np.interp(samples['time'], exact.TIMES, exact.LEVELS)),
        L2_amplitude_relative=float(np.sqrt(l2/l2norm)),
        H1_bulk_amplitude_relative=float(np.sqrt(h1/h1norm)))


def write_csv(path, rows):
    if rows:
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def save_history(path, snapshots):
    first = snapshots[0]
    for snapshot in snapshots[1:]:
        np.testing.assert_allclose(snapshot['points'], first['points'], atol=1e-12)
        np.testing.assert_allclose(snapshot['weights'], first['weights'], atol=1e-12)
    data = {key: first[key] for key in ('points', 'weights')}
    for key in ('displacement', 'gradient', 'stress', 'equivalent_plastic_strain'):
        data[key] = np.stack([s[key] for s in snapshots])
    data['times'] = np.array([s['time'] for s in snapshots])
    np.savez_compressed(path, **data)


def load_parameters(folder, args, spans):
    parameters = base.load_parameters(folder, args, spans)
    if abs(parameters['problem_data']['start_time'].GetDouble()) > 1e-14:
        raise ValueError('Cyclic loading must start from the virgin state at time zero.')
    parameters['problem_data']['end_time'].SetDouble(exact.END_TIME)
    dt = parameters['solver_settings']['time_stepping']['time_step'].GetDouble()
    count = .001/dt if dt > 0 else 0
    if count < 1 or not np.isclose(count, round(count), rtol=0, atol=1e-9):
        raise ValueError('time_step must divide 0.001 exactly, so every load reversal is reached.')
    filename = parameters['solver_settings']['material_import_settings']['materials_filename'].GetString()
    path = (folder/filename).resolve()
    material = KM.Parameters(path.read_text())
    for prop in material['properties'].values():
        law = prop['Material']['constitutive_law']['name'].GetString()
        if law != 'SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw':
            raise ValueError('Cyclic reference requires SmallStrainJ2PlasticityPlaneStrain2DThreeComponentsLaw.')
        variables = prop['Material']['Variables']
        required = dict(YOUNG_MODULUS=exact.E, POISSON_RATIO=exact.NU,
            YIELD_STRESS=exact.YIELD, ISOTROPIC_HARDENING_MODULUS=exact.H,
            EXPONENTIAL_SATURATION_YIELD_STRESS=exact.YIELD, HARDENING_EXPONENT=0.)
        for name, value in required.items():
            if not variables.Has(name) or not np.isclose(variables[name].GetDouble(), value):
                raise ValueError(f'Cyclic analytical reference requires {name}={value:g}; check {path}.')
    return parameters, material


def solve(folder, args, spans, output, capture_history=False):
    output.mkdir(parents=True, exist_ok=True)
    parameters, material = load_parameters(folder, args, spans)
    # Freeze the actual material for this run alongside the parameters.
    material_path = (output/'materials_used.json').resolve()
    material_path.write_text(material.PrettyPrintJsonString()+'\n')
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(material_path))
    (output/'parameters_used.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    metrics, snapshots = [], []

    class Analysis(base.Analysis):
        def EvaluateSolution(self, points, time):
            return exact.evaluate(points, time)

        def BoundaryTraction(self, points, time):
            return exact.boundary_traction(points, time)

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            row = errors(self.samples)
            info = self._GetSolver().GetComputingModelPart().ProcessInfo
            row['step'] = int(info[KM.STEP])
            row['iterations'] = int(info[KM.NL_ITERATION_NUMBER])
            metrics.append(row)
            self.accepted_samples = self.samples
            write_csv(output/'errors_by_time.csv', metrics)
            if capture_history:
                snapshots.append({key: value.copy() if isinstance(value, np.ndarray) else value
                                  for key, value in self.samples.items()})
            print('CYCLIC_STEP', json.dumps(row), flush=True)

    # Remove only stale products of this particular run before solving again.
    for filename in ('history.npz', 'solution.npz', 'errors_by_time.csv', 'status.json', 'contours.png'):
        (output/filename).unlink(missing_ok=True)
    with base.shared.working_directory(folder):
        Path('txt_files').mkdir(exist_ok=True)
        analysis = Analysis(KM.Model(), parameters)
        try:
            analysis.Run()
            status = dict(status='converged', time=float(analysis.time))
        except RuntimeError as error:
            status = dict(status='failed', time=float(analysis.time), error=str(error))
        counts = Counter(base.condition_name(c) for c in analysis.model['IgaModelPart'].Conditions)
    status.update(completed_steps=len(metrics), spans=list(spans), degree=args.degree)
    (output/'status.json').write_text(json.dumps(status, indent=2)+'\n')
    (output/'conditions.json').write_text(json.dumps(dict(counts), indent=2)+'\n')
    if metrics:
        np.savez_compressed(output/'solution.npz', **analysis.accepted_samples)
    if snapshots:
        save_history(output/'history.npz', snapshots)
    print('CYCLIC_RESULT', json.dumps(status), flush=True)
    return metrics, status


def reference_lines(ax, h, curves, degree):
    reference_h = np.sort(np.unique(h))
    if len(reference_h) == 1:
        reference_h = reference_h[0]*np.array([.5, 1.])
    for values, order, style in ((curves[1], degree, '--'), (curves[0], degree+1, ':')):
        valid = np.isfinite(values) & (values > 0)
        amplitude = np.exp(np.mean(np.log(values[valid])-order*np.log(h[valid]))) if valid.any() else 1.
        ax.loglog(reference_h, amplitude*reference_h**order, style, color='black', label=rf'$h^{{{order}}}$')


def plot_convergence(rows, histories, degree, output, method):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    h = np.array([r['h'] for r in rows])
    for ax, suffix, title in zip(axes, ('max', 'final'), ('Maximum over all load steps', 'Final residual state')):
        curves = [np.array([r[f'{key}_{suffix}'] for r in rows]) for key in ('L2', 'H1_bulk')]
        for values, label in zip(curves, ('L2', 'H1 bulk')):
            ax.loglog(h, values, 'o-', label=label)
        reference_lines(ax, h, curves, degree)
        ax.set(xlabel='h', ylabel='Error / first-peak reference norm', title=title)
        ax.grid(True, which='both', alpha=.3); ax.legend()
    fig.suptitle(f'{method}: cyclic annulus, p={degree}')
    fig.savefig(output/'convergence.png', dpi=180)
    detailed, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    for ax, time in zip(axes.flat, exact.TIMES[1:]):
        selected = [next(r for r in history if np.isclose(r['time'], time, atol=1e-12, rtol=0))
                    for history in histories]
        curves = [np.array([r[key] for r in selected]) for key in
                  ('L2_amplitude_relative', 'H1_bulk_amplitude_relative')]
        for values, label in zip(curves, ('L2', 'H1 bulk')):
            ax.loglog(h, values, 'o-', label=label)
        reference_lines(ax, h, curves, degree)
        ax.set(title=f't={time:g}, C/C_final={np.interp(time, exact.TIMES, exact.LEVELS):g}', xlabel='h')
        ax.grid(True, which='both', alpha=.3); ax.legend(fontsize=8)
    detailed.savefig(output/'convergence_by_time.png', dpi=180)
    return fig, detailed


def main(method, action):
    parser = argparse.ArgumentParser(description=f'Cyclic annular torsion: {method}, {action}')
    parser.add_argument('--boundary', choices=('mixed', 'dirichlet'), default='mixed')
    parser.add_argument('--parameters', help='Alternative cyclic JSON, relative to the method folder')
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--spans', nargs=2, type=int, default=[20, 20])
    # parser.add_argument('--refinements', nargs='+', type=int, default=[10, 20, 40, 80])
    parser.add_argument('--refinements', nargs='+', type=int, default=[10, 20, 40])#, 80, 160, 320])
    parser.add_argument('--dt', type=float)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--strategy', choices=('newton_raphson', 'line_search'))
    parser.add_argument('--echo-level', type=int, default=0)
    parser.add_argument('--output', type=Path, help='Override the results directory')
    parser.add_argument('--no-show', action='store_true', help='Save output without opening windows')
    parser.add_argument('--post-only', type=Path, help='Open an existing history.npz without solving')
    args = parser.parse_args()
    if min(args.degree, args.threads, *args.spans, *args.refinements) < 1 or args.echo_level < 0:
        parser.error('Degree, threads and mesh sizes must be positive; echo level must be nonnegative.')
    if args.post_only and action != 'run_and_post':
        parser.error('--post-only is available only for run_and_post_cyclic.py.')
    folder = ROOT/method
    args.parameters = args.parameters or f'ProjectParameters_{args.boundary}_cyclic.json'
    output = (args.output or folder/'results'/'cyclic'/args.boundary/f'p{args.degree}').resolve()
    KM.ParallelUtilities.SetNumThreads(args.threads)
    if action == 'run_and_post':
        from cyclic_contours import ContourViewer
        target = output/f'spans_{args.spans[0]}x{args.spans[1]}'
        if args.post_only:
            history = args.post_only.resolve()
            target = history.parent
        else:
            try:
                metrics, status = solve(folder, args, args.spans, target, capture_history=True)
            except ValueError as error:
                parser.error(str(error))
            if not metrics:
                raise RuntimeError(f'No converged load steps; see {target/"status.json"}')
            history = target/'history.npz'
            if status['status'] == 'failed':
                print(f'Solve failed: viewer contains only {len(metrics)} converged steps. See {target/"status.json"}.', flush=True)
        viewer = ContourViewer(history)
        viewer.figure.savefig(target/'contours.png', dpi=180)
        if not args.no_show:
            plt.show()
        plt.close(viewer.figure)
        return
    output.mkdir(parents=True, exist_ok=True)
    for filename in ('convergence.csv', 'convergence.png', 'convergence_by_time.png'):
        (output/filename).unlink(missing_ok=True)
    rows, histories, failures = [], [], []
    (output/'failures.json').write_text('[]\n')
    for n in sorted(set(args.refinements)):
        try:
            metrics, status = solve(folder, args, [n, n], output/f'spans_{n}x{n}')
        except ValueError as error:
            parser.error(str(error))
        if status['status'] != 'converged':
            failures.append(status)
            (output/'failures.json').write_text(json.dumps(failures, indent=2)+'\n')
            print(f'Skipping failed mesh {n}x{n}.', flush=True)
            continue
        row = dict(spans=n, h=4.4/n)
        for label, key in (('L2', 'L2_amplitude_relative'), ('H1_bulk', 'H1_bulk_amplitude_relative')):
            for suffix, value in (('max', max(m[key] for m in metrics)), ('final', metrics[-1][key])):
                name = f'{label}_{suffix}'
                row[name] = value
                row[name+'_slope'] = (float(np.log(rows[-1][name]/value)/np.log(rows[-1]['h']/row['h']))
                    if rows and min(rows[-1][name], value) > 0 else float('nan'))
        rows.append(row); histories.append(metrics)
        write_csv(output/'convergence.csv', rows)
    if rows:
        figures = plot_convergence(rows, histories, args.degree, output, method)
        if not args.no_show:
            plt.show()
        for figure in figures:
            plt.close(figure)
    print(f'Cyclic convergence: {len(rows)} completed meshes, {len(failures)} failed; {output}', flush=True)
