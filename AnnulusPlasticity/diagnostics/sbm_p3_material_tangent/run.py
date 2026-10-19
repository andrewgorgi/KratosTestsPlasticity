"""Frozen-history audits and actual Kratos utility tangent solver experiments."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
KRATOS=ROOT.parents[2]/'Kratos'
SOURCE=KRATOS/'applications/IgaApplication/custom_conditions/sbm_load_solid_condition.cpp'
LIB=KRATOS/'bin/Release/libs/libKratosIgaCore.so'
CLA=KRATOS/'bin/Release/libs/libKratosConstitutiveLawsCore.so'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    previous=HERE.parent/'sbm_p3_directional_warm_start'
    meta=json.loads((previous/'build.json').read_text())
    assert sha(SOURCE)==meta['sha256'] and sha(LIB)==meta['installed_library_sha256']
    meta['constitutive_library']=dict(path=str(CLA),sha256=sha(CLA))
    subprocess.run([sys.executable,str(HERE/'build_bridge.py')],check=True)
    (HERE/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
    for variant in ('local','utility_central','utility_forward'):
        out=HERE/variant;out.mkdir(exist_ok=True)
        for name in ('parameters_input.json','materials_input.json'):(out/name).write_bytes((previous/name).read_bytes())
        (out/'build.json').write_text(json.dumps(meta,indent=2)+'\n')
    def execute(task):
        variant,n=task;case=HERE/variant/'new'/f'spans_{n}x{n}';case.mkdir(parents=True,exist_ok=True)
        env=dict(os.environ,LD_PRELOAD=str(LIB),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/sbm-material-mpl')
        with (case/'run.log').open('w') as f:
            proc=subprocess.run([sys.executable,str(HERE/'worker.py'),'--spans',str(n),'--variant',variant],env=env,stdout=f,stderr=subprocess.STDOUT)
        if proc.returncode:raise RuntimeError(f'{variant} mesh{n}: {case}/run.log')
        result=json.loads((case/'result.json').read_text());print(variant,n,json.dumps(result['status']),flush=True);return result
    tasks=[('local',20),('local',40),('utility_central',20),('utility_central',40),('utility_forward',20)]
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(execute,tasks))
    assert sha(SOURCE)==meta['sha256'] and sha(LIB)==meta['installed_library_sha256'] and sha(CLA)==meta['constitutive_library']['sha256']
    (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
