"""Compare installed SBM Neumann condition with unit normal-dot assembly factor.
Compile an isolated real C++ condition/library; no source/install mutation.
"""
from concurrent.futures import ThreadPoolExecutor
import difflib,hashlib,json,os,re,shlex,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
KRATOS=ROOT.parents[2]/'Kratos'
SOURCE=KRATOS/'applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp'
LIB=KRATOS/'bin/Release/libs/libKratosIgaCore.so'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    original=SOURCE.read_bytes();text=original.decode();before=sha(LIB)
    # Replace precisely all three multiplications in LHS/RHS; keep normals and tractions.
    modified,count=re.subn(r'\* mTrueDotSurrogateNormal(?=[ ;])','* 1.0',text)
    assert count==3,count
    assert not re.search(r'^\s*#define ANALYTICAL_TANGENT',text,re.M)
    assert re.search(r'^\s*#define SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW',text,re.M)
    (HERE/'source_original.cpp').write_bytes(original)
    (HERE/'source_unit_factor.cpp').write_bytes(modified.encode())
    patch=''.join(difflib.unified_diff(text.splitlines(True),modified.splitlines(True),fromfile='a/applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp',tofile='b/applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp'))
    (HERE/'unit_normal_factor.patch').write_text(patch)
    temp=Path('/tmp/sbm_p3_unit_normal_factor');temp.mkdir(exist_ok=True)
    cpp=temp/'condition.cpp';obj=temp/'condition.o';library=temp/'libKratosIgaCore.so'
    cpp.write_text(modified)
    entries=json.loads((KRATOS/'build/Release/compile_commands.json').read_text())
    entry=next(e for e in entries if e['file']==str(SOURCE))
    command=shlex.split(entry['command']);command[command.index('-o')+1]=str(obj);command[-1]=str(cpp)
    with (HERE/'build.log').open('w') as log:
        subprocess.run(command,cwd=entry['directory'],stdout=log,stderr=subprocess.STDOUT,check=True)
        linkdir=KRATOS/'build/Release/applications/IgaApplication'
        command=shlex.split((linkdir/'CMakeFiles/KratosIgaCore.dir/link.txt').read_text())
        command[command.index('-o')+1]=str(library)
        command=[str(obj) if s.endswith('/sbm_load_solid_condition.cpp.o') else s for s in command]
        subprocess.run(command,cwd=linkdir,stdout=log,stderr=subprocess.STDOUT,check=True)
    meta=dict(source=str(SOURCE),sha256=hashlib.sha256(original).hexdigest(),installed_library_sha256=before,
        variants={'old':'unchanged normal dot product','new':'unit factor in all three assembly terms'},
        libraries={'old':dict(path=str(LIB),sha256=before),'new':dict(path=str(library),sha256=sha(library))})
    (HERE/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
    for name,path in [('parameters_input.json',ROOT/'sbm/ProjectParameters_mixed_cyclic.json'),('materials_input.json',ROOT/'materials.json')]:
        (HERE/name).write_bytes(path.read_bytes())
    worker=HERE.parent/'sbm_cyclic_old_vs_new/worker.py'
    def execute(task):
        mode,n=task;case=HERE/mode/f'spans_{n}x{n}';case.mkdir(parents=True,exist_ok=True)
        env=dict(os.environ,LD_PRELOAD=meta['libraries'][mode]['path'],OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/sbm-unit-mpl')
        with (case/'run.log').open('w') as log:
            process=subprocess.run([sys.executable,str(worker),'--mode',mode,'--degree','3','--spans',str(n),'--output',str(HERE)],env=env,stdout=log,stderr=subprocess.STDOUT)
        if process.returncode:raise RuntimeError(f'{mode} {n}: see {case}/run.log')
        result=json.loads((case/'result.json').read_text());print(mode,n,json.dumps(result['status']),flush=True);return result
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(execute,[(mode,n) for n in (20,40) for mode in ('old','new')]))
    (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    assert SOURCE.read_bytes()==original and sha(LIB)==before
if __name__=='__main__':main()
