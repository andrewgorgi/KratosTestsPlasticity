"""Isolated p3 sensitivity runs using the existing projection-point condition."""
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
VARIANTS = {
    'free_dt125': (.000125, -1.),
    'free_dt0625': (.0000625, -1.),
    'free_dt025': (.000025, -1.),
    'penalty_1E': (.00025, 67000.),
    'penalty_10E': (.00025, 670000.),
    'penalty_100E': (.00025, 6700000.),
    'penalty_100E_dt0625': (.0000625, 6700000.),
    'penalty_01E': (.00025, 6700.),
    'penalty_001E': (.00025, 670.),
    'penalty_1E_dt125': (.000125, 67000.),
    'penalty_1E_dt0625': (.0000625, 67000.),
}


def execute(name, spans, metadata):
    output = HERE/name
    output.mkdir(parents=True, exist_ok=True)
    dt, penalty = VARIANTS[name]
    parameters = KM.Parameters((BASE/'parameters_input.json').read_text())
    parameters['solver_settings']['time_stepping']['time_step'].SetDouble(dt)
    material = KM.Parameters((BASE/'materials_input.json').read_text())
    for prop in material['properties'].values():
        prop['Material']['Variables']['PENALTY_FACTOR'].SetDouble(penalty)
    for filename, data in [('parameters_input.json', parameters), ('materials_input.json', material)]:
        (output/filename).write_text(data.PrettyPrintJsonString()+'\n')
    (output/'build.json').write_text(json.dumps(metadata, indent=2)+'\n')
    case = output/'new'/f'spans_{spans}x{spans}'
    case.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg',
               LD_PRELOAD=metadata['libraries']['new']['path'])
    with (case/'run.log').open('w') as stream:
        process = subprocess.run([sys.executable, str(WORKER), '--mode', 'new', '--degree', '3',
                                  '--spans', str(spans), '--output', str(output)],
                                 env=env, stdout=stream, stderr=subprocess.STDOUT)
    result_path = case/'result.json'
    if process.returncode:
        result = dict(process_exit_code=process.returncode)
    else:
        result = json.loads(result_path.read_text())
    row = dict(variant=name, dt=dt, penalty=penalty, spans=spans, result=result)
    print('SENSITIVITY_CASE', json.dumps(row), flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variants', nargs='+', choices=tuple(VARIANTS), default=list(VARIANTS))
    parser.add_argument('--spans', type=int, default=20)
    args = parser.parse_args()
    metadata = json.loads((BASE/'build.json').read_text())
    source = Path(metadata['source'])
    library = source.parents[3]/'bin/Release/libs/libKratosIgaCore.so'
    for path, checksum in [(source, metadata['sha256']),
                           (library, metadata['installed_library_sha256']),
                           (Path(metadata['libraries']['new']['path']), metadata['libraries']['new']['sha256'])]:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum, path
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(lambda name: execute(name, args.spans, metadata), args.variants))
    (HERE/f'batch_{args.spans}_{"_".join(args.variants)}.json').write_text(json.dumps(rows, indent=2)+'\n')
    assert hashlib.sha256(source.read_bytes()).hexdigest() == metadata['sha256']
    assert hashlib.sha256(library.read_bytes()).hexdigest() == metadata['installed_library_sha256']
    if any('process_exit_code' in row['result'] for row in rows):
        raise RuntimeError('Worker process failed; inspect case run.log.')


if __name__ == '__main__':
    main()
