"""Check geometry, integration and solver causes with the installed p3 condition."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import KratosMultiphysics as KM

HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'sbm_neumann_current'
VARIANTS = ('diagnostic', 'newton', 'limited_ls', 'lambda05', 'lambda1', 'quadrature8', 'residual_ls', 'derivatives10')


def execute(task, metadata):
    name, n = task
    out = HERE/name
    out.mkdir(parents=True, exist_ok=True)
    parameters = KM.Parameters((BASE/'parameters_input.json').read_text())
    if name=='newton':
        parameters['solver_settings']['solving_strategy_settings']['type'].SetString('newton_raphson')
    for modeler in parameters['modelers'].values():
        settings = modeler['Parameters']
        if name in ('lambda05', 'lambda1') and settings.Has('lambda_inner'):
            settings['lambda_inner'].SetDouble(.5 if name=='lambda05' else 1.)
        if name=='quadrature8' and settings.Has('element_condition_list'):
            entry = settings['element_condition_list'][0]
            entry.AddInt('number_of_integration_points_per_span', 8)
        if name=='derivatives10' and settings.Has('element_condition_list'):
            for entry in settings['element_condition_list'].values():
                if entry['type'].GetString()=='condition':
                    entry['shape_function_derivatives_order'].SetInt(10)
    (out/'parameters_input.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    (out/'materials_input.json').write_bytes((BASE/'materials_input.json').read_bytes())
    (out/'build.json').write_text(json.dumps(metadata, indent=2)+'\n')
    case = out/'new'/f'spans_{n}x{n}'
    case.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg',
               LD_PRELOAD=metadata['libraries']['new']['path'])
    with (case/'run.log').open('w') as log:
        process = subprocess.run([sys.executable, str(HERE/'worker.py'), '--variant', name,
                                 '--spans', str(n)], stdout=log, stderr=subprocess.STDOUT, env=env)
    if process.returncode:
        raise RuntimeError(f'Worker failed: {name}, mesh {n}; see {case}/run.log')
    result = json.loads((case/'result.json').read_text())
    print('GLOBAL_CAUSE_CASE', name, n, json.dumps(result['status']), flush=True)
    return dict(variant=name, **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variants', nargs='+', choices=VARIANTS, default=list(VARIANTS))
    parser.add_argument('--spans', type=int, nargs='+', default=[20])
    args = parser.parse_args()
    metadata = json.loads((BASE/'build.json').read_text())
    for path, sha in [(Path(metadata['source']), metadata['sha256']),
                      (Path(metadata['libraries']['new']['path']), metadata['libraries']['new']['sha256'])]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==sha, path
    tasks = [(name, n) for n in args.spans for name in args.variants]
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda task: execute(task, metadata), tasks))
    (HERE/f'batch_{"_".join(args.variants)}_{"_".join(map(str,args.spans))}.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__=='__main__':
    main()
