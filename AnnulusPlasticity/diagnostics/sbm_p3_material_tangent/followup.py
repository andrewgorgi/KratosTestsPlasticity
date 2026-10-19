"""Check cycle continuation and perturbation sensitivity after the recovered p3 step."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'build.json').read_text());lib=meta['libraries']['new'];assert hashlib.sha256(Path(lib['path']).read_bytes()).hexdigest()==lib['sha256']
def execute(task):
    tag,n,extra=task;out=HERE/tag;out.mkdir(exist_ok=True)
    for name in ('parameters_input.json','materials_input.json','build.json'):(out/name).write_bytes((HERE/'utility_central'/name).read_bytes())
    case=out/'new'/f'spans_{n}x{n}';case.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,LD_PRELOAD=lib['path'],OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/sbm-material-mpl')
    with (case/'run.log').open('w') as f:
        p=subprocess.run([sys.executable,str(HERE/'worker.py'),'--spans',str(n),'--variant','utility_central','--tag',tag,*extra],env=env,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'{tag}: {case}/run.log')
    result=json.loads((case/'result.json').read_text());print(tag,n,json.dumps(result['status']),flush=True);return result
with ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(execute,[('utility_central_full',20,['--full']),('utility_central_h1e10',20,['--h','1e-10'])]))
(HERE/'followup_results.json').write_text(json.dumps(results,indent=2)+'\n')
