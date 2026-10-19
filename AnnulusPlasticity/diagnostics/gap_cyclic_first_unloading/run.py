import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import cyclic_study as study

parser = argparse.ArgumentParser()
parser.add_argument('--strategy', choices=['newton_raphson', 'line_search'], required=True)
parser.add_argument('--quadrature', type=int, choices=[0,4], required=True)
parser.add_argument('--audit-jacobian', action='store_true')
args = parser.parse_args()
study.KM.ParallelUtilities.SetNumThreads(1)
out = Path(__file__).parent/f'{args.strategy}_q{args.quadrature or "default"}'
out.mkdir(exist_ok=True)
source = ROOT/'gap_sbm'
parameters = study.KM.Parameters((Path(__file__).parent/'parameters_input.json').read_text())
parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(Path(__file__).parent/'materials_input.json'))
for modeler in parameters['modelers'].values():
    p = modeler['Parameters']
    if p.Has('input_filename'):
        p['input_filename'].SetString(str((source/p['input_filename'].GetString()).resolve()))
parameters['solver_settings']['solving_strategy_settings']['type'].SetString(args.strategy)
entry = parameters['modelers'][3]['Parameters']['element_condition_list'][0]
if args.quadrature:
    entry.AddInt('number_of_integration_points_per_span', args.quadrature)
(out/'input.json').write_text(parameters.PrettyPrintJsonString())
# Only the first load ramp and first unloading step are needed for this audit.
study.exact.END_TIME = .00125
original = study.base.Analysis.InitializeSolutionStep
audit = []
def initialize(self):
    original(self)
    part = self.model['IgaModelPart']
    time = float(part.ProcessInfo[study.KM.TIME])
    c = study.exact.C_FINAL*np.interp(time, study.exact.TIMES, study.exact.LEVELS)
    force_error = 0.
    loads = 0
    for condition in part.Conditions:
        if 'GapSbmLoadSolidCondition' not in study.base.condition_name(condition):
            continue
        center = condition.GetGeometry().Center()
        r = np.hypot(center.X, center.Y)
        expected = c*np.array([center.Y,-center.X])/r**3
        actual = np.array(condition.GetValue(study.KM.FORCE))[:2]
        force_error = max(force_error,float(np.max(abs(actual-expected))))
        loads += 1
    body_error = max(np.linalg.norm(np.array(e.GetValue(study.KM.BODY_FORCE))) for e in part.Elements)
    audit.append(dict(time=time, load_conditions=loads, max_traction_error=force_error, body_force=body_error))
study.base.Analysis.InitializeSolutionStep = initialize
if args.audit_jacobian:
    from audit_global_cycle import audit as audit_jacobian
    original_initialize = study.base.Analysis.Initialize
    def initialize_analysis(self):
        original_initialize(self)
        solver = self._GetSolver()
        checked_solve = solver.SolveSolutionStep
        def checked_with_audit():
            try:
                return checked_solve()
            except RuntimeError:
                destination = out/'jacobian'
                destination.mkdir(exist_ok=True)
                audit_jacobian(solver, destination)
                from entity_audit import audit_entities
                audit_entities(solver, destination)
                raise
        solver.SolveSolutionStep = checked_with_audit
    study.base.Analysis.Initialize = initialize_analysis
settings = SimpleNamespace(parameters='input.json', boundary='mixed', degree=2, dt=None, echo_level=1)
metrics, status = study.solve(out, settings, [10,10], out/'results')
(out/'boundary_audit.json').write_text(json.dumps(audit,indent=2))
print('CHECK_RESULT', json.dumps(dict(status=status, boundary_audit=audit)),flush=True)
