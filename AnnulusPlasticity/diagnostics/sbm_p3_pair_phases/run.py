"""Reproduce the native p3 failures and audit surrogate/true phase pairs."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import KratosMultiphysics as KM

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def execute(n, metadata):
    case = HERE/'new'/f'spans_{n}x{n}'
    case.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg',
               LD_PRELOAD=metadata['libraries']['new']['path'])
    with (case/'run.log').open('w') as log:
        process = subprocess.run([sys.executable, str(HERE/'worker.py'), '--spans', str(n)],
                                 env=env, stdout=log, stderr=subprocess.STDOUT)
    if process.returncode:
        raise RuntimeError(f'Worker failed, mesh {n}; see {case}/run.log')
    data = json.loads((case/'result.json').read_text())
    print('PAIR_PHASE_CASE', n, json.dumps(data['status']), flush=True)
    return data


def main():
    metadata = json.loads((HERE.parent/'sbm_neumann_current'/'build.json').read_text())
    for path, expected in [(Path(metadata['source']), metadata['sha256']),
                           (Path(metadata['libraries']['new']['path']), metadata['installed_library_sha256'])]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, path
    (HERE/'build.json').write_text(json.dumps(metadata, indent=2)+'\n')
    for filename, source in [('parameters_input.json', ROOT/'sbm'/'ProjectParameters_mixed_cyclic.json'),
                             ('materials_input.json', ROOT/'materials.json')]:
        (HERE/filename).write_text(KM.Parameters(source.read_text()).PrettyPrintJsonString()+'\n')
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda n:execute(n, metadata), (20, 40)))
    (HERE/'results.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__=='__main__': main()
