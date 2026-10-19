"""One isolated SBM mesh, using a frozen input and one existing macro branch."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import cyclic_study as study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('old', 'new'), required=True)
    parser.add_argument('--spans', type=int, required=True)
    parser.add_argument('--degree', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--boundary', choices=('mixed', 'dirichlet'), default='mixed')
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args()
    output = args.output.resolve()
    metadata = json.loads((output/'build.json').read_text())
    library = metadata['libraries'][args.mode]
    assert os.environ.get('LD_PRELOAD') == library['path']
    assert hashlib.sha256(Path(library['path']).read_bytes()).hexdigest() == library['sha256']
    study.KM.ParallelUtilities.SetNumThreads(1)
    case = output/args.mode/f'spans_{args.spans}x{args.spans}'
    case.mkdir(parents=True, exist_ok=True)
    parameters = study.KM.Parameters((output/'parameters_input.json').read_text())
    for modeler in parameters['modelers'].values():
        settings = modeler['Parameters']
        if settings.Has('input_filename'):
            path = (ROOT/'sbm'/settings['input_filename'].GetString()).resolve()
            settings['input_filename'].SetString(str(path))
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(
        str(output/'materials_input.json'))
    (case/'input.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    settings = SimpleNamespace(parameters='input.json', boundary=args.boundary, degree=args.degree,
                               dt=None, echo_level=0)
    start = time.monotonic()
    metrics, status = study.solve(case, settings, [args.spans]*2, case/'results')
    result = dict(mode=args.mode, degree=args.degree, spans=args.spans, seconds=time.monotonic()-start,
                  status=status, metrics=metrics)
    (case/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print('COMPARISON_CASE', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
