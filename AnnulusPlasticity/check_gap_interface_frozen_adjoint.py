"""Check the frozen interface adjoint with the installed Kratos build.

Runs in an isolated directory, preserving production files. Audits the actual
global and entity Jacobians before committing selected cyclic load steps.
"""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

import cyclic_study as study
from audit_global_cycle import audit

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'diagnostics/gap_cyclic_first_unloading'))
from entity_audit import audit_entities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spans', type=int, default=10)
    parser.add_argument('--strategy', choices=('newton_raphson', 'line_search'), default='newton_raphson')
    parser.add_argument('--end-time', type=float, default=study.exact.END_TIME)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'diagnostics/gap_interface_frozen_adjoint/check')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    study.KM.ParallelUtilities.SetNumThreads(1)
    # Reuse the exact pre-fix input to make the before/after check comparable.
    fixture = ROOT/'diagnostics/gap_cyclic_first_unloading'
    parameters = study.KM.Parameters((fixture/'parameters_input.json').read_text())
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(
        str(fixture/'materials_input.json'))
    parameters['solver_settings']['solving_strategy_settings']['type'].SetString(args.strategy)
    (out/'input.json').write_text(parameters.PrettyPrintJsonString()+'\n')

    load_parameters = study.load_parameters
    def load(folder, settings, spans):
        parameters, material = load_parameters(folder, settings, spans)
        parameters['problem_data']['end_time'].SetDouble(args.end_time)
        return parameters, material
    study.load_parameters = load

    initialize = study.base.Analysis.Initialize
    checks = []
    selected = (.0005, .00125, .003, .006)
    def initialize_and_instrument(self):
        initialize(self)
        solver = self._GetSolver()
        solve_checked = solver.SolveSolutionStep
        def solve_and_check():
            failure = None
            try:
                converged = solve_checked()
            except RuntimeError as error:
                failure = error
            time = float(solver.GetComputingModelPart().ProcessInfo[study.KM.TIME])
            if failure or any(np.isclose(time, t, atol=1e-12, rtol=0) for t in selected):
                destination = out/f'jacobian_t{time:.5f}'
                destination.mkdir(exist_ok=True)
                audit(solver, destination)
                audit_entities(solver, destination)
                global_data = json.loads((destination/'audit.json').read_text())
                entity_data = json.loads((destination/'entity_audit.json').read_text())
                row = dict(time=time,
                    global_fd_error=max(item['error'] for item in global_data['jacobian']),
                    matrix_reconstruction_error=entity_data['matrix_reconstruction_error'],
                    repeat_error=entity_data['repeat_error'], by_type=entity_data['by_type'])
                checks.append(row)
                (out/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
                print('FROZEN_ADJOINT_CHECK', json.dumps(row), flush=True)
            if failure:
                raise failure
            return converged
        solver.SolveSolutionStep = solve_and_check
    study.base.Analysis.Initialize = initialize_and_instrument
    settings = SimpleNamespace(parameters='input.json', boundary='mixed', degree=2, dt=None, echo_level=0)
    metrics, status = study.solve(out, settings, [args.spans]*2, out/'results')
    (out/'result.json').write_text(json.dumps(dict(status=status, checks=checks), indent=2)+'\n')
    assert status['status'] == 'converged', status
    for check in checks:
        assert check['matrix_reconstruction_error'] < 1e-10, check
        assert check['repeat_error'] < 1e-10, check
        assert check['global_fd_error'] < 1e-5, check


if __name__ == '__main__':
    main()
