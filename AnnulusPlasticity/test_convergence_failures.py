"""Convergence continues past failed meshes and plots only valid results."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import annulus_study as study


class ConvergenceFailureTests(unittest.TestCase):
    def test_failed_mesh_is_skipped_and_successful_meshes_are_plotted(self):
        def solve(folder, args, spans):
            if spans[0] == 10:
                raise RuntimeError('Newton did not converge')
            return {'mesh': spans[0]}, study.KM.Parameters('{}'), {'SbmLoadSolidCondition': 8}

        def errors(samples):
            return dict(L2_relative=1./samples['mesh']**3,
                        H1_bulk_relative=1./samples['mesh']**2)

        with TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(study, 'ROOT', root), \
                 patch.object(study, 'solve', side_effect=solve) as mocked_solve, \
                 patch.object(study, 'errors', side_effect=errors), \
                 patch.object(study.plt, 'show') as show, \
                 patch('sys.argv', ['convergence.py', '--refinements', '9', '10', '11']):
                study.main('sbm', 'convergence')
            output = root/'sbm/results/dirichlet/p2'
            self.assertEqual(mocked_solve.call_count, 3)
            with (output/'convergence.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([int(r['spans_x']) for r in rows], [9, 11])
            self.assertAlmostEqual(float(rows[1]['L2_relative_slope']), 3.)
            failures = json.loads((output/'failures.json').read_text())
            self.assertEqual(len(failures), 1)
            self.assertEqual(failures[0]['spans_x'], 10)
            self.assertIn('Newton did not converge', failures[0]['error'])
            self.assertTrue((output/'convergence.png').exists())
            show.assert_called_once()

    def test_default_range_and_all_failed_run_clear_stale_summary(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            output = root/'sbm/results/dirichlet/p2'
            output.mkdir(parents=True)
            for filename in ('convergence.csv', 'convergence.png'):
                (output/filename).write_text('Previous run')
            with patch.object(study, 'ROOT', root), \
                 patch.object(study, 'solve', side_effect=RuntimeError('Failed')) as solve, \
                 patch.object(study.plt, 'show') as show, \
                 patch('sys.argv', ['convergence.py']):
                study.main('sbm', 'convergence')
            self.assertEqual([call.args[2] for call in solve.call_args_list],
                             [(n, n) for n in range(9, 60)])
            self.assertEqual(len(json.loads((output/'failures.json').read_text())), 51)
            self.assertFalse((output/'convergence.csv').exists())
            self.assertFalse((output/'convergence.png').exists())
            show.assert_not_called()

    def test_single_run_still_reports_failure(self):
        with TemporaryDirectory() as directory:
            with patch.object(study, 'ROOT', Path(directory)), \
                 patch.object(study, 'solve', side_effect=RuntimeError('Failed')), \
                 patch('sys.argv', ['run_and_post.py']):
                with self.assertRaisesRegex(RuntimeError, 'Failed'):
                    study.main('sbm', 'run_and_post')


if __name__ == '__main__':
    unittest.main()
