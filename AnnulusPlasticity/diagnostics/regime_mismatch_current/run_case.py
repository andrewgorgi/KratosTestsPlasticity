"""Run one p2 mesh with the separately instrumented IGA library preloaded.

Read frozen input files next to this script. Negative NL_ITERATION_NUMBER in
REGIME_PAIR records denotes a final trial snapshot before commit, or failure.
"""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
import annulus_study as study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True)
    parser.add_argument('--strategy', choices=('line_search', 'newton_raphson'))
    args = parser.parse_args()
    study.KM.ParallelUtilities.SetNumThreads(1)
    case = OUT/(args.strategy or 'current')/f'n{args.mesh}'
    case.mkdir(parents=True, exist_ok=True)
    records = []

    class Analysis(study.Analysis):
        def Initialize(self):
            super().Initialize()
            solver = self._GetSolver()
            solve = solver.SolveSolutionStep

            def snapshot(converged):
                part = solver.GetComputingModelPart()
                info = part.ProcessInfo
                iteration = int(info[study.KM.NL_ITERATION_NUMBER])
                record = dict(time=float(info[study.KM.TIME]),
                              iteration=iteration, converged=converged)
                try:
                    info[study.KM.NL_ITERATION_NUMBER] = -max(iteration, 1)
                    for condition in part.Conditions:
                        if 'SbmLoadSolidCondition' in study.condition_name(condition):
                            matrix = study.KM.Matrix()
                            condition.CalculateLeftHandSide(matrix, info)
                except Exception as error:
                    record['snapshot_error'] = str(error)
                finally:
                    info[study.KM.NL_ITERATION_NUMBER] = iteration
                records.append(record)
                (case/'steps.json').write_text(json.dumps(records, indent=2)+'\n')

            def solve_and_snapshot():
                try:
                    converged = solve()
                except RuntimeError:
                    snapshot(False)
                    raise
                snapshot(bool(converged))
                return converged

            solver.SolveSolutionStep = solve_and_snapshot

    settings = SimpleNamespace(parameters=str(OUT/'parameters_input.json'),
        boundary='mixed', degree=2, dt=None, strategy=args.strategy, echo_level=0)
    parameters = study.load_parameters(ROOT/'sbm', settings, [args.mesh]*2)
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(
        str(OUT/'materials_input.json'))
    (case/'parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    with study.shared.working_directory(ROOT/'sbm'):
        analysis = Analysis(study.KM.Model(), parameters)
        try:
            analysis.Run()
            result = dict(status='converged', time=float(analysis.time))
        except Exception as error:
            result = dict(status='failed', time=float(analysis.time), error=str(error))
    result['spans'] = args.mesh
    result['strategy'] = parameters['solver_settings']['solving_strategy_settings']['type'].GetString()
    result['steps'] = records
    (case/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print('REGIME_RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
