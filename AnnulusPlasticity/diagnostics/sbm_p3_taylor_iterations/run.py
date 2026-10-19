"""Freeze current inputs/library; compare native p3 solves with/without replay audit."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,subprocess,sys
from pathlib import Path
import KratosMultiphysics as KM
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=Path('/home/agorgi/Documents/GECKO/Kratos/applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp')
LIBRARY=Path('/home/agorgi/Documents/GECKO/Kratos/bin/Release/libs/libKratosIgaCore.so')
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run_case(task):
    mode,n=task;out=HERE/mode
    case=out/'new'/f'spans_{n}x{n}';case.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/sbm-taylor-mpl',LD_PRELOAD=str(LIBRARY))
    with (case/'run.log').open('w') as f:
        proc=subprocess.run([sys.executable,str(HERE/'worker.py'),'--spans',str(n),'--mode',mode],env=env,stdout=f,stderr=subprocess.STDOUT)
    if proc.returncode:raise RuntimeError(f'{mode} n={n}: see {case}/run.log')
    result=json.loads((case/'result.json').read_text())
    print(mode,n,json.dumps(result['status']),flush=True)
    return result

def main():
    source=SOURCE.read_text()
    assert '#define SBM_USE_PROJECTION_POINT_CONSTITUTIVE_LAW' in source
    assert 'std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder)' in source
    assert '// #define ANALYTICAL_TANGENT' in source
    meta=dict(source=str(SOURCE),sha256=digest(SOURCE),installed_library_sha256=digest(LIBRARY),libraries={'new':dict(path=str(LIBRARY),sha256=digest(LIBRARY))},strain_taylor_order=2)
    (HERE/'source_snapshot.cpp').write_bytes(SOURCE.read_bytes())
    (HERE/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
    for mode in ('baseline','audit'):
        out=HERE/mode;out.mkdir(exist_ok=True)
        (out/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
        for filename,src in [('parameters_input.json',ROOT/'sbm/ProjectParameters_mixed_cyclic.json'),('materials_input.json',ROOT/'materials.json')]:
            (out/filename).write_text(KM.Parameters(src.read_text()).PrettyPrintJsonString()+'\n')
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(run_case,[(mode,n) for n in (20,40) for mode in ('baseline','audit')]))
    assert digest(SOURCE)==meta['sha256'] and digest(LIBRARY)==meta['installed_library_sha256']
    (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
