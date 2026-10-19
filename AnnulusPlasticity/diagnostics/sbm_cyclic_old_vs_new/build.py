"""Build the two existing macro branches in /tmp; never install either library."""
import argparse
import hashlib
import json
import re
import shlex
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
KRATOS = ROOT.parents[2] / 'Kratos'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE)
    parser.add_argument('--inputs', type=Path, help='Reuse frozen parameters/materials from an earlier comparison')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    build = Path('/tmp')/output.name
    source = KRATOS/'applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp'
    original = source.read_bytes()
    text = original.decode()
    assert not re.search(r'^\s*#define ANALYTICAL_TANGENT', text, re.M)
    entries = json.loads((KRATOS/'build/Release/compile_commands.json').read_text())
    entry = next(e for e in entries if e['file'] == str(source))
    linkdir = KRATOS/'build/Release/applications/IgaApplication'
    metadata = {'source': str(source), 'sha256': hashlib.sha256(original).hexdigest(),
                'libraries': {}, 'installed_library_sha256': hashlib.sha256((KRATOS/'bin/Release/libs/libKratosIgaCore.so').read_bytes()).hexdigest()}
    for mode in ('old', 'new'):
        folder = build/mode
        folder.mkdir(parents=True, exist_ok=True)
        content, count = re.subn(
            r'^(?:\s*//\s*)?#define SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW\s*$',
            ('#define ' if mode == 'new' else '// #define ')
            + 'SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW', text, flags=re.M)
        assert count == 1
        cpp, obj, library = folder/'condition.cpp', folder/'condition.o', folder/'libKratosIgaCore.so'
        cpp.write_text(content)
        command = shlex.split(entry['command'])
        command[command.index('-o')+1] = str(obj)
        command[-1] = str(cpp)
        with (output/f'build_{mode}.log').open('w') as log:
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
    for name, path in [('parameters_input.json', ROOT/'sbm/ProjectParameters_mixed_cyclic.json'),
                       ('materials_input.json', ROOT/'materials.json')]:
        input_path = args.inputs/name if args.inputs else path
        (output/name).write_bytes(input_path.read_bytes())
    (output/'build.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
