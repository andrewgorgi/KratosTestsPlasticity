"""Summarize and plot existing cyclic audit records without rerunning Kratos."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    local, global_results = [], []
    all_rows = []
    for mode in ('on', 'off'):
        rows = json.loads((HERE/f'local_{mode}.json').read_text())
        all_rows.extend(rows)
        for law in ('linear', 'saturation'):
            for degree in (1, 2):
                for spans in sorted({r['spans'] for r in rows}):
                    selected = [r for r in rows if r['law'] == law
                                and r['degree'] == degree and r['spans'] == spans]
                    keys = [k for k in selected[0] if k.endswith('_error')]
                    local.append(dict(mode=mode, law=law, degree=degree, spans=spans,
                                      states=len(selected), **{k: max(r[k] for r in selected) for k in keys}))
    for path in sorted((HERE/'global').glob('*/*/result.json')):
        record = json.loads(path.read_text())
        global_results.append(record)
    summary = dict(local=local, global_results=global_results,
                   material_control=json.loads((HERE/'material_control.json').read_text()),
                   reference_validation=json.loads((HERE/'reference_validation.json').read_text()))
    (HERE/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for mode in ('on', 'off'):
        for spans, color in ((20, 'tab:blue'), (40, 'tab:orange')):
            selected = [r for r in all_rows if r['mode'] == mode and r['law'] == 'linear'
                        and r['degree'] == 2 and r['spans'] == spans]
            label = f'macro {mode}, n={spans}'
            for ax, key in zip(axes[0], ('stress_history_relative_error', 'finite_difference_error')):
                values = [max(r[key] for r in selected if r['step'] == step) for step in range(1, 9)]
                ax.semilogy(range(1, 9), np.maximum(values, 1e-16),
                            '--' if mode == 'off' else '-', marker='o', color=color, label=label)
    axes[0, 0].set(title='Local p2: projected stress vs independent history',
                   xlabel='Load step', ylabel='Maximum relative stress error')
    axes[0, 1].set(title='Local p2: assembled tangent vs finite differences',
                   xlabel='Load step', ylabel='Maximum directional Jacobian error')
    for record in global_results:
        if record['degree'] != 2:
            continue
        steps = record['steps']
        color = 'tab:blue' if record['spans'] == 20 else 'tab:orange'
        axes[1, 0].plot([r['time'] for r in steps], [r['L2_amplitude_relative'] for r in steps],
            '--' if record['mode'] == 'off' else '-', color=color, marker='.',
            label=f"macro {record['mode']}, n={record['spans']}")
        if record['status'] == 'failed':
            axes[1, 0].plot(steps[-1]['time'], steps[-1]['L2_amplitude_relative'],
                            'x', color=color, ms=10)
    axes[1, 0].set(title='Global p2: cyclic displacement error (x: last step before failure)',
                   xlabel='Time', ylabel='L2 error / reference at first peak')
    from reference import LEVELS, TIMES
    axes[1, 1].plot(TIMES, LEVELS, 'ko-')
    axes[1, 1].set(title='Prescribed torque history', xlabel='Time', ylabel='C(t) / C_final')
    for ax in axes.flat:
        ax.grid(True, alpha=.3)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.legend(fontsize=8)
    fig.savefig(HERE/'comparison.png', dpi=180)
    plt.close(fig)
    print('Written', HERE/'summary.json', HERE/'comparison.png')


if __name__ == '__main__':
    main()
