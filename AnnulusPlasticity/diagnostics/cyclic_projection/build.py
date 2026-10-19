"""Build the two existing macro branches in /tmp; never install either library."""
import hashlib
import json
import re
import shlex
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
KRATOS = ROOT.parents[2] / 'Kratos'
BUILD = Path('/tmp/sbm_cyclic_projection')


def main():
    source = KRATOS/'applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp'
    original = source.read_bytes()
    text = original.decode()
    assert not re.search(r'^\s*#define ANALYTICAL_TANGENT', text, re.M)
    entries = json.loads((KRATOS/'build/Release/compile_commands.json').read_text())
    entry = next(e for e in entries if e['file'] == str(source))
    linkdir = KRATOS/'build/Release/applications/IgaApplication'
    metadata = {'source': str(source), 'sha256': hashlib.sha256(original).hexdigest(),
                'libraries': {}}
    for mode in ('off', 'on'):
        folder = BUILD/mode
        folder.mkdir(parents=True, exist_ok=True)
        content, count = re.subn(
            r'^(?:\s*//\s*)?#define SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW\s*$',
            ('#define ' if mode == 'on' else '// #define ')
            + 'SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW', text, flags=re.M)
        assert count == 1
        cpp, obj, library = folder/'condition.cpp', folder/'condition.o', folder/'libKratosIgaCore.so'
        cpp.write_text(content)
        command = shlex.split(entry['command'])
        command[command.index('-o')+1] = str(obj)
        command[-1] = str(cpp)
        with (HERE/f'build_{mode}.log').open('w') as log:
            subprocess.run(command, cwd=entry['directory'], stdout=log,
                           stderr=subprocess.STDOUT, check=True)
            command = shlex.split((linkdir/'CMakeFiles/KratosIgaCore.dir/link.txt').read_text())
            command[command.index('-o')+1] = str(library)
            command = [str(obj) if s.endswith('/sbm_load_solid_condition.cpp.o') else s
                       for s in command]
            subprocess.run(command, cwd=linkdir, stdout=log, stderr=subprocess.STDOUT, check=True)
        metadata['libraries'][mode] = {'path': str(library),
            'sha256': hashlib.sha256(library.read_bytes()).hexdigest()}
        print('Built', mode, library, flush=True)
    assert source.read_bytes() == original
    for name, path in [('parameters_input.json', ROOT/'sbm/ProjectParameters_mixed.json'),
                       ('materials_input.json', ROOT/'materials.json')]:
        (HERE/name).write_bytes(path.read_bytes())
    (HERE/'build.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
