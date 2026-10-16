"""Simulazione breve del caso reale; interrompe se Newton non converge."""
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import study

degree = int(sys.argv[1]) if len(sys.argv)>1 else 2
threads = int(sys.argv[2]) if len(sys.argv)>2 else 4
study.KM.ParallelUtilities.SetNumThreads(threads)
iterations = []


class Checked(study.Analysis):
    def Initialize(self):
        super().Initialize()
        solver = self._GetSolver()
        native = solver.SolveSolutionStep

        def solve_checked():
            success = native()
            iterations.append(dict(time=self.time,converged=bool(success),
                                   iterations=solver.GetComputingModelPart().ProcessInfo[study.KM.NL_ITERATION_NUMBER]))
            if not success:
                raise RuntimeError('Newton non convergente')
            return success

        solver.SolveSolutionStep = solve_checked


parameters = study.load_parameters(ROOT/'gap_sbm',SimpleNamespace(
    parameters=None,field='cosine_boundary',physics='plastic',spans=[8,4],degree=degree))
result = dict(degree=degree,threads=threads,spans=[8,4])
with tempfile.TemporaryDirectory(prefix='cut_analytical_run_',dir='/tmp') as temporary:
    root = Path(temporary)
    (root/'geometry').symlink_to(ROOT/'geometry',target_is_directory=True)
    (root/'materials_plastic.json').symlink_to(ROOT/'materials_plastic.json')
    folder=root/'gap_sbm';folder.mkdir();(folder/'txt_files').mkdir()
    with study.working_directory(folder):
        try:
            analysis = Checked(study.KM.Model(),parameters)
            analysis.Run()
            result.update(success=True,**study.errors(analysis.samples,'cosine_boundary','plastic'))
        except Exception as error:
            result.update(success=False,error=str(error))
result['iterations'] = iterations
out = Path(__file__).resolve().parent
(out/f'gap_analytical_p{degree}_threads{threads}_results.json').write_text(json.dumps(result,indent=2)+'\n')
(out/f'gap_analytical_p{degree}_threads{threads}_parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
print('GAP_ANALYTICAL_RUN',json.dumps(result),flush=True)
if not result['success']:
    raise SystemExit(1)
