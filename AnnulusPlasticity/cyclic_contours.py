"""Matplotlib contour viewer with synchronized load-step and field controls."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.path import Path as PlotPath
from matplotlib.patches import PathPatch
import matplotlib.tri as mtri
from matplotlib.widgets import RadioButtons, Slider
import numpy as np

import cyclic_analytic as exact


class ContourViewer:
    fields = ('u_x', 'u_y', '|u|', 'sigma_xx', 'sigma_yy', 'sigma_xy',
              'equivalent_plastic_strain')

    def __init__(self, filename):
        self.filename = Path(filename)
        with np.load(filename) as archive:
            self.data = {key: archive[key] for key in archive.files}
        self.times = self.data['times']
        self.field = 'equivalent_plastic_strain'
        self.reference = [exact.evaluate(self.data['points'], t) for t in self.times]
        points = self.data['points']
        self.tri = mtri.Triangulation(*points.T)
        centers = points[self.tri.triangles].mean(axis=1)
        r = np.linalg.norm(centers, axis=1)
        self.tri.set_mask((r < exact.RI) | (r > exact.RO))
        self.figure = plt.figure(figsize=(15, 8))
        self.axes = [self.figure.add_axes([left, .39, .25, .51]) for left in (.045, .375, .705)]
        self.colorbars, self.artists = [], []
        self.clips = []
        for ax in self.axes:
            bar = self.figure.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap='viridis'),
                                      ax=ax, fraction=.046, pad=.035)
            self.colorbars.append(bar)
            ax.set(xlabel='x', ylabel='y', aspect='equal', xlim=(-2.05, 2.05), ylim=(-2.05, 2.05))
            angle = np.linspace(0, 2*np.pi, 401)
            for radius in (exact.RI, exact.RO):
                ax.plot(radius*np.cos(angle), radius*np.sin(angle), 'k-', lw=.7)
            # Cut contours at the exact physical circles, including triangles
            # whose centroids lie in the annulus but whose edges cross the hole.
            vertices = np.concatenate([np.column_stack((radius*np.cos(a), radius*np.sin(a)))
                                       for radius, a in ((exact.RO, angle), (exact.RI, angle[::-1]))])
            codes = np.full(len(vertices), PlotPath.LINETO, dtype=np.uint8)
            codes[0] = codes[len(angle)] = PlotPath.MOVETO
            codes[len(angle)-1] = codes[-1] = PlotPath.CLOSEPOLY
            self.clips.append(PathPatch(PlotPath(vertices, codes), transform=ax.transData))
        self.title = self.figure.suptitle('')
        self.radio = RadioButtons(self.figure.add_axes([.02, .035, .21, .27]),
                                  self.fields, active=len(self.fields)-1)
        for label in self.radio.labels:
            label.set_fontsize(9)
        load_ax = self.figure.add_axes([.31, .16, .62, .14])
        load_ax.plot(exact.TIMES, exact.LEVELS, 'k.-')
        load_ax.set(xlabel='Tempo', ylabel='C / C_final')
        load_ax.grid(True, alpha=.3)
        self.load_marker, = load_ax.plot([], [], 'ro')
        self.load_line = load_ax.axvline(self.times[0], color='tab:red', alpha=.4)
        self.slider = Slider(self.figure.add_axes([.31, .06, .62, .035]),
                             'Step', 1, max(2, len(self.times)), valinit=1, valstep=1, valfmt='%0.0f')
        if len(self.times) == 1:
            self.slider.set_active(False)
        self.slider.on_changed(self.update)
        self.radio.on_clicked(self.select_field)
        self.update(1)

    def values(self, field):
        u = self.data['displacement']
        s = self.data['stress']
        if field == 'equivalent_plastic_strain':
            numerical = self.data[field]
            analytical = np.stack([r[3] for r in self.reference])
        elif field.startswith('sigma_'):
            column = ('sigma_xx', 'sigma_yy', 'sigma_xy').index(field)
            numerical = s[:, :, column]
            analytical = np.stack([r[2][:, column] for r in self.reference])
        elif field == '|u|':
            numerical = np.linalg.norm(u, axis=2)
            analytical = np.stack([np.linalg.norm(r[0], axis=1) for r in self.reference])
        else:
            column = ('u_x', 'u_y').index(field)
            numerical = u[:, :, column]
            analytical = np.stack([r[0][:, column] for r in self.reference])
        return numerical, analytical

    def select_field(self, field):
        self.field = field
        self.update(self.slider.val)

    def update(self, step):
        index = min(max(int(step)-1, 0), len(self.times)-1)
        numerical, analytical = self.values(self.field)
        difference = numerical-analytical
        # Fixed ranges across the whole saved history make frames comparable.
        low = min(numerical.min(), analytical.min())
        high = max(numerical.max(), analytical.max())
        if high-low < 1e-14:
            high = low+1e-12
        error_scale = max(float(np.max(abs(difference))), 1e-12)
        for artist in self.artists:
            if hasattr(artist, 'remove'):
                artist.remove()
            else:
                for collection in artist.collections:
                    collection.remove()
        self.artists = []
        for ax, colorbar, clip, values, label, limits, cmap in zip(
                self.axes, self.colorbars, self.clips, (numerical[index], analytical[index], difference[index]),
                ('Numerico', 'Analitico', 'Errore numerico - analitico'),
                ((low, high), (low, high), (-error_scale, error_scale)),
                ('viridis', 'viridis', 'coolwarm')):
            norm = Normalize(*limits)
            artist = ax.tricontourf(self.tri, values, levels=np.linspace(*limits, 31),
                                    cmap=cmap, norm=norm, extend='both')
            if hasattr(artist, 'set_clip_path'):
                artist.set_clip_path(clip)
            else:
                for collection in artist.collections:
                    collection.set_clip_path(clip)
            self.artists.append(artist)
            colorbar.update_normal(ScalarMappable(norm=norm, cmap=cmap))
            ax.set_title(label)
        time = self.times[index]
        level = float(np.interp(time, exact.TIMES, exact.LEVELS))
        self.load_marker.set_data([time], [level])
        self.load_line.set_xdata([time, time])
        partial = ' — ciclo parziale' if self.times[-1] < exact.END_TIME-1e-12 else ''
        self.title.set_text(f'{self.field} — step {index+1}/{len(self.times)}, t={time:g}, C/C_final={level:g}{partial}')
        self.figure.canvas.draw_idle()
