"""Compare Neumann quadrature points and the analytical plastic front, mesh 40."""
import csv
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[1]))
from analytic import evaluate, FINAL_PLASTIC_RADIUS, END_TIME

time = .001
front = FINAL_PLASTIC_RADIUS*np.sqrt(time/END_TIME)
angle = np.linspace(0, 2*np.pi, 600)
fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
for ax, value in zip(axes, (.5, 0.)):
    points = {}
    with (OUT/f'lambda_{value:g}'/'current/n40/pairs.csv').open() as stream:
        for row in csv.DictReader(stream):
            points.setdefault(row['condition'], [float(row['x_s']), float(row['y_s'])])
    points = np.array(list(points.values()))
    plastic = evaluate(points, time)[3] > 0
    ax.plot(np.cos(angle), np.sin(angle), 'k-', label='True boundary, r=1')
    ax.plot(front*np.cos(angle), front*np.sin(angle), '--', color='tab:orange',
            label=f'Analytical plastic front, r={front:.4f}')
    for mask, color, label in [(plastic, 'tab:red', 'Surrogate: plastic reference'),
                               (~plastic, 'tab:blue', 'Surrogate: elastic reference')]:
        if mask.any():
            ax.scatter(points[mask, 0], points[mask, 1], s=7, color=color, label=label)
    ax.set(title=f'lambda_inner={value:g}, n=40, t={time:g}',
           xlabel='x', ylabel='y', xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
    ax.set_aspect('equal')
    ax.grid(True, alpha=.2)
    ax.legend(fontsize=8, loc='upper right')
fig.savefig(OUT/'surrogate_geometry.png', dpi=180)
plt.close(fig)
