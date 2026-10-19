"""Replace inner Neumann by exact Dirichlet, preserving the mixed-case geometry."""
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
BASE = HERE.parent/'sbm_cyclic_old_vs_new_p3'
WORKER = HERE.parent/'sbm_cyclic_old_vs_new'/'worker.py'
VARIANTS = {'free': (-1., .00025), 'penalty_100E': (6700000., .00025),
            'free_dt125': (-1., .000125)}


def execute(task, metadata):
    name, spans = task
    penalty, dt = VARIANTS[name]
    output = HERE/name
    case = output/'new'/f'spans_{spans}x{spans}'
    case.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg',
               LD_PRELOAD=metadata['libraries']['new']['path'])
    if not (case/'result.json').exists():
        with (case/'run.log').open('w') as stream:
            process = subprocess.run([sys.executable, str(WORKER), '--mode', 'new', '--degree', '3',
                '--boundary', 'dirichlet', '--spans', str(spans), '--output', str(output)],
                env=env, stdout=stream, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(f'{name} {spans}: worker exit {process.returncode}; see {case}/run.log')
    data = json.loads((case/'result.json').read_text())
    print('DIRICHLET_CASE', name, spans, json.dumps(data['status']), flush=True)
    return dict(variant=name, penalty=penalty, dt=dt, **data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spans', type=int, nargs='+', default=[20, 40])
    parser.add_argument('--variants', nargs='+', choices=tuple(VARIANTS), default=['free', 'penalty_100E'])
    args = parser.parse_args()
    metadata = json.loads((BASE/'build.json').read_text())
    source = Path(metadata['source'])
    library = source.parents[3]/'bin/Release/libs/libKratosIgaCore.so'
    for path, checksum in [(source, metadata['sha256']),
                           (library, metadata['installed_library_sha256']),
                           (Path(metadata['libraries']['new']['path']), metadata['libraries']['new']['sha256'])]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==checksum, path
    for name in args.variants:
        output = HERE/name
        output.mkdir(parents=True, exist_ok=True)
        parameters = KM.Parameters((BASE/'parameters_input.json').read_text())
        changed = 0
        for modeler in parameters['modelers'].values():
            settings = modeler['Parameters']
            if settings.Has('link_layer_to_condition_name'):
                for layer in settings['link_layer_to_condition_name'].values():
                    if layer['condition_name'].GetString()=='SbmLoadSolidCondition':
                        layer['condition_name'].SetString('SbmSolidCondition')
                        changed += 1
        assert changed==1
        penalty, dt = VARIANTS[name]
        parameters['solver_settings']['time_stepping']['time_step'].SetDouble(dt)
        material = KM.Parameters((BASE/'materials_input.json').read_text())
        for prop in material['properties'].values():
            prop['Material']['Variables']['PENALTY_FACTOR'].SetDouble(penalty)
        for filename, data in [('parameters_input.json', parameters), ('materials_input.json', material)]:
            (output/filename).write_text(data.PrettyPrintJsonString()+'\n')
        (output/'build.json').write_text(json.dumps(metadata, indent=2)+'\n')
    tasks = [(name, n) for n in args.spans for name in args.variants]
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda task:execute(task, metadata), tasks))
    (HERE/f'batch_{"_".join(args.variants)}_{"_".join(map(str,args.spans))}.json').write_text(
        json.dumps(results, indent=2)+'\n')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==metadata['sha256']
    assert hashlib.sha256(library.read_bytes()).hexdigest()==metadata['installed_library_sha256']


if __name__=='__main__':
    main()
