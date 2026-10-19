"""Partial cycles contribute only at times actually solved to convergence."""
import csv
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import postprocess


class PartialCycleComparisonTests(unittest.TestCase):
    def test_p1_reference_orders(self):
        rows = [dict(mode='new', h=h, L2_amplitude_relative=h**2,
                     H1_bulk_amplitude_relative=h) for h in (.2, .1)]
        figure, ax = postprocess.plt.subplots()
        postprocess.references(ax, rows, degree=1)
        self.assertEqual(set(ax.get_legend_handles_labels()[1]), {r'$h^{2}$', r'$h^{1}$'})
        postprocess.plt.close(figure)

    def test_partial_states_and_slopes_use_equal_times(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'execution.json').write_text(json.dumps(dict(refinements=[10, 20])))
            for mode in ('old', 'new'):
                for n in (10, 20):
                    times = [.001, .00125, .006] if mode == 'new' else ([.001] if n == 10 else [.001, .00125])
                    h = 4.4/n
                    factor = 2 if mode == 'old' else 1
                    metrics = [dict(time=t, L2_amplitude_relative=factor*h**3,
                                    H1_bulk_amplitude_relative=factor*h**2) for t in times]
                    destination = root/mode/f'spans_{n}x{n}'
                    destination.mkdir(parents=True)
                    (destination/'result.json').write_text(json.dumps(dict(mode=mode, spans=n, seconds=1,
                        status=dict(status='converged' if mode == 'new' else 'failed', time=.006 if mode == 'new' else times[-1]+.00025),
                        metrics=metrics)))
            with patch.object(postprocess, 'HERE', root):
                postprocess.main()
            with (root/'comparison.csv').open() as stream:
                comparison = list(csv.DictReader(stream))
            self.assertEqual([(int(r['spans']), float(r['time'])) for r in comparison],
                             [(10, .001), (20, .001), (20, .00125)])
            self.assertTrue(all(float(r['L2_amplitude_relative_old_over_new']) == 2. for r in comparison))
            with (root/'slopes_by_time.csv').open() as stream:
                slopes = list(csv.DictReader(stream))
            first_peak = next(r for r in slopes if r['mode'] == 'old' and float(r['time']) == .001 and int(r['spans']) == 20)
            self.assertAlmostEqual(float(first_peak['L2_amplitude_relative_slope']), 3.)
            unloading = next(r for r in slopes if r['mode'] == 'old' and float(r['time']) == .00125)
            self.assertTrue(math.isnan(float(unloading['L2_amplitude_relative_slope'])))
            summary = json.loads((root/'summary.json').read_text())
            self.assertEqual(summary['common_converged_states'], 3)
            self.assertEqual(summary['error_states'], 9)


if __name__ == '__main__':
    unittest.main()
