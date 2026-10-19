"""Cyclic reference, residual normalization, skipped meshes and slider behavior."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import cyclic_analytic as exact
import cyclic_study as study
from cyclic_contours import ContourViewer


class CyclicTests(unittest.TestCase):
    def samples(self, time):
        angle = np.linspace(0, 2*np.pi, 20, endpoint=False)
        r, theta = np.meshgrid(np.linspace(1.02, 1.98, 5), angle)
        points = np.column_stack(((r*np.cos(theta)).ravel(), (r*np.sin(theta)).ravel()))
        u, gradient, stress, alpha = exact.evaluate(points, time)
        return dict(points=points, weights=np.ones(len(points)), displacement=u,
                    gradient=gradient, stress=stress, equivalent_plastic_strain=alpha, time=time)

    def test_reversal_reference_and_residual_normalization(self):
        self.assertLess(exact.verify()['scalar_history_error'], 1e-13)
        for time in (.001, .002, .003, .004, .005, .006):
            samples = self.samples(time)
            self.assertEqual(study.errors(samples)['L2_amplitude_relative'], 0.)
            self.assertEqual(study.errors(samples)['H1_bulk_amplitude_relative'], 0.)
        errors = []
        for time in (.001, .002, .006):
            samples = self.samples(time)
            samples['displacement'] += .001
            errors.append(study.errors(samples)['L2_amplitude_relative'])
        np.testing.assert_allclose(errors, errors[0], rtol=1e-13)
        residual = self.samples(.002)
        self.assertGreater(np.linalg.norm(residual['displacement']), 0.)
        np.testing.assert_allclose(residual['stress'], 0., atol=1e-12)

    def test_slider_changes_step_and_field_with_fixed_color_scales(self):
        with TemporaryDirectory() as temporary:
            filename = Path(temporary)/'history.npz'
            study.save_history(filename, [self.samples(t) for t in (.001, .002, .003)])
            viewer = ContourViewer(filename)
            initial_limits = viewer.artists[0].get_clim()
            viewer.slider.set_val(3)
            self.assertIn('t=0.003', viewer.title.get_text())
            self.assertEqual(viewer.artists[0].get_clim(), initial_limits)
            viewer.radio.set_active(0)
            self.assertEqual(viewer.field, 'u_x')
            self.assertEqual(len(viewer.artists), 3)
            self.assertEqual(len(viewer.figure.axes), 9)
            plt.close(viewer.figure)

    def test_convergence_skips_failed_cycle_and_keeps_both_reference_orders(self):
        def solve(folder, args, spans, output):
            if spans[0] == 12:
                return [], dict(status='failed', time=.00125, error='Newton failed')
            metrics = [dict(time=float(t), L2_amplitude_relative=1/spans[0]**3,
                            H1_bulk_amplitude_relative=1/spans[0]**2) for t in exact.TIMES[1:]]
            return metrics, dict(status='converged', time=exact.END_TIME)

        with TemporaryDirectory() as temporary:
            output = Path(temporary)
            with patch.object(study, 'solve', side_effect=solve), \
                 patch.object(study.plt, 'show') as show, \
                 patch('sys.argv', ['convergence_cyclic.py', '--output', str(output),
                                    '--refinements', '10', '12', '20']):
                study.main('sbm', 'convergence')
            with (output/'convergence.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([int(r['spans']) for r in rows], [10, 20])
            self.assertAlmostEqual(float(rows[1]['L2_max_slope']), 3.)
            self.assertAlmostEqual(float(rows[1]['H1_bulk_max_slope']), 2.)
            self.assertEqual(len(json.loads((output/'failures.json').read_text())), 1)
            self.assertTrue((output/'convergence_by_time.png').exists())
            show.assert_called_once()

    def test_time_step_must_reach_turning_points(self):
        from types import SimpleNamespace
        args = SimpleNamespace(parameters='ProjectParameters_mixed_cyclic.json',
                               boundary='mixed', degree=2, dt=.0003)
        with self.assertRaisesRegex(ValueError, 'divide 0.001'):
            study.load_parameters(study.ROOT/'sbm', args, [10, 10])


if __name__ == '__main__':
    unittest.main()
