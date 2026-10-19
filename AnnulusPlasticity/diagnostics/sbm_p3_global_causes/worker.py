"""Existing-condition experiments; optional global-state and residual-merit audit."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker as original_worker


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant', required=True)
    parser.add_argument('--spans', type=int, required=True)
    args = parser.parse_args()
    output = HERE/args.variant
    if args.variant=='limited_ls':
        from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
        KM = original_worker.study.KM
        def configured(self):
            settings = KM.Parameters(json.dumps(dict(max_iteration=self.settings['max_iteration'].GetInt(),
                compute_reactions=False, reform_dofs_at_each_step=False, move_mesh_flag=False,
                first_alpha_value=.5, second_alpha_value=1., min_alpha=1e-4, max_alpha=1.,
                max_line_search_iterations=10, line_search_tolerance=.5)))
            strategy = KM.LineSearchStrategy(self.GetComputingModelPart(), self._GetScheme(),
                self._GetConvergenceCriterion(), self._GetBuilderAndSolver(), settings)
            strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings['use_old_stiffness_in_first_iteration'].GetBool())
            (output/'line_search_effective_settings.json').write_text(settings.PrettyPrintJsonString()+'\n')
            return strategy
        MechanicalSolver._create_line_search_strategy = configured
    if args.variant=='diagnostic':
        from global_audit import install
        install(original_worker.study, output/f'audit_n{args.spans}')
    if args.variant=='residual_ls':
        from residual_solver import install
        install(original_worker.study, output/f'newton_trace_n{args.spans}')
    sys.argv = [__file__, '--mode', 'new', '--degree', '3', '--spans', str(args.spans), '--output', str(output)]
    original_worker.main()


if __name__=='__main__':
    main()
