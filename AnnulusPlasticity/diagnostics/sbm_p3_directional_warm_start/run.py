"""Run p3 directional-Jacobian and L2-projected warm-start controls."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
KRATOS=ROOT.parents[2]/'Kratos'
SOURCE=KRATOS/'applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp'
LIB=KRATOS/'bin/Release/libs/libKratosIgaCore.so'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    meta=dict(source=str(SOURCE),sha256=sha(SOURCE),installed_library_sha256=sha(LIB),libraries={'new':dict(path=str(LIB),sha256=sha(LIB))})
    (HERE/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
    (HERE/'source_snapshot.cpp').write_bytes(SOURCE.read_bytes())
    for name,path in [('parameters_input.json',ROOT/'sbm/ProjectParameters_mixed_cyclic.json'),('materials_input.json',ROOT/'materials.json')]:
        (HERE/name).write_bytes(path.read_bytes())
    def execute(n):
        case=HERE/'new'/f'spans_{n}x{n}';case.mkdir(parents=True,exist_ok=True)
        env=dict(os.environ,LD_PRELOAD=str(LIB),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/sbm-directional-mpl')
        with (case/'run.log').open('w') as f:
            process=subprocess.run([sys.executable,str(HERE/'worker.py'),'--spans',str(n)],env=env,stdout=f,stderr=subprocess.STDOUT)
        if process.returncode:raise RuntimeError(f'Mesh {n}: see {case}/run.log')
        result=json.loads((case/'result.json').read_text());print(n,json.dumps(result['status']),flush=True);return result
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(execute,(20,40)))
    assert sha(SOURCE)==meta['sha256'] and sha(LIB)==meta['installed_library_sha256']
    (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
