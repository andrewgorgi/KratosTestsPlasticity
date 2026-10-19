"""Run one p2 mesh with the separately instrumented IGA library preloaded.

Read frozen input files next to this script. Negative NL_ITERATION_NUMBER in
REGIME_PAIR records denotes a final trial snapshot before commit, or failure.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from types import SimpleNamespace

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
import annulus_study as study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True)
    parser.add_argument('--lambda-inner', type=float, default=0.)
    parser.add_argument('--strategy', choices=('line_search', 'newton_raphson'))
    args = parser.parse_args()
    study.KM.ParallelUtilities.SetNumThreads(1)
    case = OUT/f'lambda_{args.lambda_inner:g}'/(args.strategy or 'current')/f'n{args.mesh}'
    case.mkdir(parents=True, exist_ok=True)
    records = []
    geometry_info = {}

    class Analysis(study.Analysis):
        def Initialize(self):
            super().Initialize()
            part = self.model['IgaModelPart']
            conditions = [c for c in part.Conditions
                          if 'SbmLoadSolidCondition' in study.condition_name(c)]
            centers = np.array([[c.GetGeometry().Center().X, c.GetGeometry().Center().Y]
                                for c in conditions])
            radii = np.linalg.norm(centers, axis=1)
            analytic_alpha = study.evaluate(centers, .001)[3]
            geometry_info.update(points=len(conditions),
                surrogate_radius_min=float(radii.min()),
                surrogate_radius_max=float(radii.max()),
                analytic_elastic_surrogate_points_at_first_step=int((analytic_alpha==0).sum()))
            (case/'geometry.json').write_text(json.dumps(geometry_info, indent=2)+'\n')
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
    study.shared.geometry_settings(parameters)['lambda_inner'].SetDouble(args.lambda_inner)
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(
        str(OUT/'materials_input.json'))
    (case/'parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    with study.shared.working_directory(ROOT/'sbm'):
        analysis = Analysis(study.KM.Model(), parameters)
        try:
            analysis.Run()
            result = dict(status='converged', **study.errors(analysis.samples))
            np.savez(case/'solution.npz', **analysis.samples)
        except Exception as error:
            result = dict(status='failed', time=float(analysis.time), error=str(error))
    result['lambda_inner'] = args.lambda_inner
    result['geometry'] = geometry_info
    result['spans'] = args.mesh
    result['strategy'] = parameters['solver_settings']['solving_strategy_settings']['type'].GetString()
    result['steps'] = records
    (case/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print('REGIME_RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
