"""Reproduce the two existing macro branches and all cyclic diagnostic tests.

Run: python3 AnnulusPlasticity/diagnostics/cyclic_projection/run.py
Only temporary libraries are built. Installed Kratos and production inputs are
not changed. Solver failures are test observations recorded in result.json.
"""
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def execute(script, arguments=(), mode=None, log=None):
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg')
    if mode:
        environment['LD_PRELOAD'] = f'/tmp/sbm_cyclic_projection/{mode}/libKratosIgaCore.so'
    command = [sys.executable, str(HERE/script), *arguments]
    print('Running', script, *arguments, flush=True)
    if log:
        with (HERE/log).open('w') as stream:
            subprocess.run(command, env=environment, stdout=stream,
                           stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(command, env=environment, check=True)


def main():
    execute('build.py')
    execute('material_control.py', log='material_control.log')
    for mode in ('on', 'off'):
        execute('local.py', ['--mode', mode], mode=mode, log=f'local_{mode}.log')
        for degree, spans in ((1, 20), (2, 20), (2, 40)):
            execute('global_case.py', ['--mode', mode, '--degree', str(degree), '--spans', str(spans)],
                    mode=mode, log=f'global_{mode}_p{degree}_n{spans}.log')
    execute('summarize.py')


if __name__ == '__main__':
    main()
