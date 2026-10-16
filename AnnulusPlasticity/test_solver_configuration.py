"""Solver overrides and diagnostic preservation, without an FE solve."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import annulus_study as study


class SolverConfigurationTests(unittest.TestCase):
    def args(self,**overrides):
        values = dict(parameters=None,boundary='dirichlet',degree=1,dt=None,
                      strategy=None,echo_level=None,gap_load_from_force=False)
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_gap_defaults_and_explicit_overrides(self):
        folder = study.ROOT/'gap_sbm'
        parameters = study.load_parameters(folder,self.args(),[10,10])
        settings = parameters['solver_settings']
        self.assertEqual(settings['solving_strategy_settings']['type'].GetString(),'line_search')
        self.assertEqual(settings['time_stepping']['time_step'].GetDouble(),.001)
        override = study.load_parameters(folder,self.args(dt=.0001,strategy='newton_raphson',echo_level=1),[20,20])
        self.assertEqual(override['solver_settings']['solving_strategy_settings']['type'].GetString(),'newton_raphson')
        self.assertEqual(override['solver_settings']['time_stepping']['time_step'].GetDouble(),.0001)
        self.assertEqual(override['solver_settings']['echo_level'].GetInt(),1)
        for key in ('max_iteration','displacement_relative_tolerance','residual_relative_tolerance'):
            self.assertEqual(settings[key].GetDouble(),override['solver_settings'][key].GetDouble())

    def test_failed_solve_preserves_settings(self):
        source = study.ROOT/'gap_sbm/ProjectParameters_dirichlet.json'
        class FailingAnalysis:
            def __init__(self,*args):
                pass
            def Run(self):
                raise RuntimeError('Nonlinear solve failed in regression test')
        with TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'ProjectParameters_dirichlet.json').write_text(source.read_text())
            with patch.object(study,'Analysis',FailingAnalysis):
                with self.assertRaisesRegex(RuntimeError,'Diagnostics:'):
                    study.solve(folder,self.args(),[10,10])
            failure = folder/'results/dirichlet/p1/failed_spans_10x10'
            info = json.loads((failure/'failure.json').read_text())
            self.assertEqual(info['strategy'],'line_search')
            self.assertEqual(info['dt'],.001)
            self.assertIn('Nonlinear solve failed',info['error'])
            self.assertTrue((failure/'parameters_used.json').exists())


if __name__ == '__main__':
    unittest.main()
