"""Compare existing SBM macro branches without changing the installed build.

Default meshes: 10,20,40,80. Larger meshes can be requested with --refinements.
Numerical failures remain observations; process/tool failures are reported apart.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def execute(mode, spans, metadata, output, degree):
    case = output/mode/f'spans_{spans}x{spans}'
    case.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg',
                       LD_PRELOAD=metadata['libraries'][mode]['path'])
    with (case/'run.log').open('w') as stream:
        result = subprocess.run([sys.executable, str(HERE/'worker.py'), '--mode', mode,
                                 '--spans', str(spans), '--degree', str(degree), '--output', str(output)], env=environment,
                                stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        return dict(mode=mode, spans=spans, process_exit_code=result.returncode)
    data = json.loads((case/'result.json').read_text())
    print('CASE_COMPLETE', mode, spans, data['status']['status'], data['status']['time'], flush=True)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refinements', type=int, nargs='+', default=[10, 20, 40, 80])
    parser.add_argument('--degree', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if min(args.refinements) < 1:
        parser.error('Refinements must be positive.')
    output = (args.output or (HERE if args.degree == 2 else HERE.parent/f'sbm_cyclic_old_vs_new_p{args.degree}')).resolve()
    output.mkdir(parents=True, exist_ok=True)
    build_command = [sys.executable, str(HERE/'build.py'), '--output', str(output)]
    if args.degree != 2:
        build_command += ['--inputs', str(HERE)]
    subprocess.run(build_command, check=True)
    metadata = json.loads((output/'build.json').read_text())
    failures = []
    for spans in sorted(set(args.refinements)):
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(execute, mode, spans, metadata, output, args.degree) for mode in ('old', 'new')]
            for future in futures:
                failure = future.result()
                if failure:
                    failures.append(failure)
        (output/'process_failures.json').write_text(json.dumps(failures, indent=2)+'\n')
    source = Path(metadata['source'])
    assert hashlib.sha256(source.read_bytes()).hexdigest() == metadata['sha256']
    library = source.parents[3]/'bin/Release/libs/libKratosIgaCore.so'
    assert hashlib.sha256(library.read_bytes()).hexdigest() == metadata['installed_library_sha256']
    (output/'execution.json').write_text(json.dumps(dict(degree=args.degree, refinements=sorted(set(args.refinements)),
        source_unchanged=True, installed_library_unchanged=True), indent=2)+'\n')
    subprocess.run([sys.executable, str(HERE/'postprocess.py'), '--output', str(output)], check=True)
    if failures:
        raise RuntimeError(f'{len(failures)} processes failed; see process_failures.json.')


if __name__ == '__main__':
    main()
