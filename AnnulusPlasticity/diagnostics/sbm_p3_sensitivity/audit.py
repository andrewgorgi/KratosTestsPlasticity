"""Check the p3 residual Jacobian without changing any condition."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'sbm_cyclic_old_vs_new_p3'
sys.path.insert(0, str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker
from audit_global_cycle import audit

output = HERE/'jacobian_check'
output.mkdir(exist_ok=True)
for name in ('build.json', 'parameters_input.json', 'materials_input.json'):
    (output/name).write_bytes((BASE/name).read_bytes())
original = worker.study.base.Analysis.Initialize


def initialize(self):
    original(self)
    solver = self._GetSolver()
    solve = solver.SolveSolutionStep

    def solve_and_audit():
        try:
            result = solve()
        except RuntimeError:
            destination = output/'failed_state'
            destination.mkdir(exist_ok=True)
            audit(solver, destination)
            raise
        if any(abs(float(self.time)-time) < 1e-12 for time in (.00025, .0005)):
            destination = output/('elastic_state' if float(self.time) < .0003 else 'converged_state')
            destination.mkdir(exist_ok=True)
            audit(solver, destination)
        return result

    solver.SolveSolutionStep = solve_and_audit


worker.study.base.Analysis.Initialize = initialize
sys.argv = [__file__, '--mode', 'new', '--degree', '3', '--spans', '20', '--output', str(output)]
worker.main()
