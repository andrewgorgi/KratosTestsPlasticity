"""Build a diagnostic Python bridge to the real Kratos tangent utility in /tmp."""
import hashlib,json,shlex,subprocess,sysconfig
from pathlib import Path
HERE=Path(__file__).resolve().parent
KRATOS=HERE.parents[1].parents[2]/'Kratos'
entries=json.loads((KRATOS/'build/Release/compile_commands.json').read_text())
e=next(x for x in entries if x['file'].endswith('constitutive_laws_python_application.cpp'))
folder=Path('/tmp/sbm_material_tangent_bridge');folder.mkdir(exist_ok=True)
obj=folder/'bridge.o';library=folder/('sbm_tangent_bridge'+sysconfig.get_config_var('EXT_SUFFIX'))
cmd=shlex.split(e['command']);cmd[cmd.index('-o')+1]=str(obj);cmd[-1]=str(HERE/'bridge.cpp')
with (HERE/'build_bridge.log').open('w') as f:
    subprocess.run(cmd,cwd=e['directory'],stdout=f,stderr=subprocess.STDOUT,check=True)
    subprocess.run(['/usr/bin/g++','-shared','-fopenmp',str(obj),'-L'+str(KRATOS/'bin/Release/libs'),'-Wl,-rpath,'+str(KRATOS/'bin/Release/libs'),'-lKratosConstitutiveLawsCore','-lKratosCore','-o',str(library)],stdout=f,stderr=subprocess.STDOUT,check=True)
(HERE/'bridge.json').write_text(json.dumps(dict(path=str(library),directory=str(folder),sha256=hashlib.sha256(library.read_bytes()).hexdigest()),indent=2)+'\n')
print(library,flush=True)
