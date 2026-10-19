"""Compare spatial errors at equal times, retaining converged partial cycles."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
KEYS = ('L2_amplitude_relative', 'H1_bulk_amplitude_relative')
LABELS = {'old': 'Vecchia', 'new': 'Nuova'}
COLORS = {'old': 'tab:orange', 'new': 'tab:blue'}


def write_csv(filename, rows):
    if not rows:
        return
    with (HERE/filename).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def references(ax, rows, degree=2):
    for key, order, style in ((KEYS[0], degree+1, ':'), (KEYS[1], degree, '--')):
        chosen = [r for r in rows if r['mode'] == 'new' and r[key] > 0]
        if not chosen:
            continue
        h = np.array([r['h'] for r in chosen])
        values = np.array([r[key] for r in chosen])
        scale = np.exp(np.mean(np.log(values)-order*np.log(h)))
        span = np.unique(np.sort(h))
        if len(span) == 1:
            span = span[0]*np.array([.5, 1.])
        ax.loglog(span, scale*span**order, style, color='black', alpha=.7, label=rf'$h^{{{order}}}$')


def convergence_panel(ax, time, rows, degree=2):
    chosen = [r for r in rows if abs(r['time']-time) < 1e-10]
    for mode in ('new', 'old'):
        samples = sorted([r for r in chosen if r['mode'] == mode], key=lambda r: r['h'], reverse=True)
        for key, label, marker in zip(KEYS, ('L2', 'H1 bulk'), ('o', 's')):
            ax.loglog([r['h'] for r in samples], [r[key] for r in samples],
                      color=COLORS[mode], marker=marker, linestyle='-' if mode == 'new' else '--',
                      label=f'{LABELS[mode]}: {label}')
        if not samples:
            ax.text(.03, .05 if mode == 'old' else .13, f'{LABELS[mode]}: nessuna mesh convergente',
                    transform=ax.transAxes, color=COLORS[mode], fontsize=8)
    references(ax, chosen, degree)
    ax.set(title=f't = {time:g}', xlabel='h', ylabel='Errore / norma al primo picco')
    ax.grid(True, which='both', alpha=.3)
    ax.legend(fontsize=7)


def main(output=None):
    global HERE
    if output is not None:
        HERE = Path(output).resolve()
    execution = json.loads((HERE/'execution.json').read_text())
    degree = execution.get('degree', 2)
    cases, states, availability = [], [], []
    by_state = {}
    for n in execution['refinements']:
        for mode in ('old', 'new'):
            filename = HERE/mode/f'spans_{n}x{n}/result.json'
            if not filename.exists():
                continue
            case = json.loads(filename.read_text())
            cases.append(case)
            status, metrics = case['status'], case['metrics']
            availability.append(dict(mode=mode, spans=n, h=4.4/n, status=status['status'],
                completed_steps=len(metrics), last_converged_time=metrics[-1]['time'] if metrics else 0.,
                failed_time=status['time'] if status['status'] == 'failed' else '',
                seconds=case['seconds']))
            for metric in metrics:
                row = dict(mode=mode, spans=n, h=4.4/n, cycle_status=status['status'], **metric)
                states.append(row)
                by_state[(mode, n, round(row['time'], 10))] = row
    write_csv('availability.csv', availability)
    write_csv('errors_by_time.csv', states)
    comparison = []
    for row in states:
        if row['mode'] != 'old':
            continue
        new = by_state.get(('new', row['spans'], round(row['time'], 10)))
        if new is None:
            continue
        pair = dict(spans=row['spans'], h=row['h'], time=row['time'])
        for key in KEYS:
            pair[key+'_old'] = row[key]
            pair[key+'_new'] = new[key]
            pair[key+'_old_over_new'] = row[key]/new[key] if new[key] > 0 else float('nan')
        comparison.append(pair)
    write_csv('comparison.csv', comparison)
    slopes = []
    for time in sorted({round(r['time'], 10) for r in states}):
        for mode in ('new', 'old'):
            samples = sorted([r for r in states if r['mode'] == mode and abs(r['time']-time) < 1e-10],
                             key=lambda r: r['h'], reverse=True)
            previous = None
            for row in samples:
                entry = dict(mode=mode, time=time, spans=row['spans'], h=row['h'])
                for key in KEYS:
                    entry[key] = row[key]
                    entry[key+'_slope'] = (float(np.log(previous[key]/row[key])/np.log(previous['h']/row['h']))
                        if previous and min(previous[key], row[key]) > 0 else float('nan'))
                slopes.append(entry)
                previous = row
    write_csv('slopes_by_time.csv', slopes)
    figure, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    for ax, time in zip(axes.flat, (.001, .00125, .002, .006)):
        convergence_panel(ax, time, states, degree)
    figure.suptitle(f'SBM p{degree}: storia al surrogate (vecchia) e al punto proiettato (nuova)')
    figure.savefig(HERE/'comparison.png', dpi=180)
    plt.close(figure)
    initial_times = sorted({round(r['time'], 10) for r in states})[:4]
    if initial_times:
        figure, axes = plt.subplots(1, len(initial_times), figsize=(6*len(initial_times), 5),
                                    squeeze=False, constrained_layout=True)
        for ax, time in zip(axes.flat, initial_times):
            convergence_panel(ax, time, states, degree)
        figure.suptitle(f'SBM p{degree}: primi passi convergenti')
        figure.savefig(HERE/'initial_convergence.png', dpi=180)
        plt.close(figure)
    figure, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
    for ax, time in zip(axes.flat, (.0005, .001, .00125, .002, .003, .004, .005, .006)):
        convergence_panel(ax, time, states, degree)
    figure.savefig(HERE/'convergence_by_time.png', dpi=180)
    plt.close(figure)
    figure, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for index, row in enumerate(availability):
        ax.barh(index, row['last_converged_time']*1000, color=COLORS[row['mode']], alpha=.75)
        if row['failed_time'] != '':
            ax.plot(row['failed_time']*1000, index, 'rx', ms=8)
        else:
            ax.plot(row['last_converged_time']*1000, index, 'ko', ms=4)
    ax.set_yticks(range(len(availability)),
                  [f"{r['spans']}x{r['spans']} — {LABELS[r['mode']]}" for r in availability])
    ax.set(xlabel='Tempo [ms]', title='Ultimo passo convergente; x rosso = passo fallito', xlim=(0, 6.3))
    ax.grid(True, axis='x', alpha=.3)
    figure.savefig(HERE/'availability.png', dpi=180)
    plt.close(figure)
    summary = dict(degree=degree, availability=availability, common_converged_states=len(comparison),
                   error_states=len(states), meshes=execution['refinements'])
    (HERE/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('COMPARISON_SUMMARY', json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE)
    main(parser.parse_args().output)
