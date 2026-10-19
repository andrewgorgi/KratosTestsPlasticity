"""Plot the constitutive-regime discordance along actual Newton iterations."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
for ax, n, time, outcome in zip(axes.flat, (10, 20, 40, 80),
                               (.001, .001, .001, .002),
                               ('converged', 'failed', 'failed', 'failed')):
    with (OUT/'newton_raphson'/f'n{n}'/'iterations.csv').open() as stream:
        rows = [r for r in csv.DictReader(stream)
                if float(r['time']) == time and r['final_snapshot'] == 'False']
    x = [int(r['iteration']) for r in rows]
    ax.plot(x, [int(r['elastic_s_plastic_t']) for r in rows], 'o-',
            markersize=2, linewidth=1, label='Elastic surrogate / plastic true')
    ax.plot(x, [int(r['plastic_s_elastic_t']) for r in rows], 'o-',
            markersize=2, linewidth=1, label='Plastic surrogate / elastic true')
    ax.set(title=f'Newton p2, {n}x{n}, t={time:g}: {outcome}',
           xlabel='Newton iteration', ylabel='Discordant Neumann pairs')
    ax.grid(True, alpha=.3)
    if n == 80:
        ax.axvline(2.5, color='gray', linestyle=':', linewidth=1)
        ax.text(.98, .9, 'Iterations 3-100: no discordance',
                transform=ax.transAxes, ha='right', va='top')
axes[0, 0].legend(fontsize=8)
fig.savefig(OUT/'regime_cycles.png', dpi=180)
plt.close(fig)
