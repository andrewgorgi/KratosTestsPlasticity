"""Repeat projected starts in independent processes, without a prior failed solve."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'build.json').read_text());lib=meta['libraries']['new']
assert hashlib.sha256(Path(lib['path']).read_bytes()).hexdigest()==lib['sha256']
out=HERE/'fresh';out.mkdir(exist_ok=True)
for name in ('build.json','parameters_input.json','materials_input.json'):(out/name).write_bytes((HERE/name).read_bytes())
def run(n):
    case=out/'new'/f'spans_{n}x{n}';case.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,LD_PRELOAD=lib['path'],OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/sbm-directional-mpl')
    with (case/'run.log').open('w') as f:
        p=subprocess.run([sys.executable,str(HERE/'worker.py'),'--spans',str(n),'--fresh'],env=env,stdout=f,stderr=subprocess.STDOUT)
    assert p.returncode==0,(n,case)
    result=json.loads((case/'result.json').read_text());print(n,json.dumps(result['status']),flush=True)
    return result
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,(20,40)))
(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
