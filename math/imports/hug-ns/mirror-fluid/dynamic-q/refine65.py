"""Four additional spatial controls authorized after the first 26-run result."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import hashlib,os,time
import study as s
JOBS=[dict(id=name+'-n65',case=name,n=65,half=False,mirror=False) for name,*_ in s.TEST]
def main():
    lock=s.HERE/'refine65.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    try:
        for path,sha in s.read(s.HERE/'training-freeze.json')['files'].items():
            assert hashlib.sha256((s.HERE/path).read_bytes()).hexdigest()==sha,path
        protocol=dict(authorization='User said keep going after results were reported; announced the four 65-cubed controls.',
          jobs=JOBS,end_time=s.END,starting_fields='Exactly the four original unseen starts, sampled at65 cubed.',
          fits='Frozen training coefficients unchanged. These controls are not fitting data.',
          source_hashes=s.hashes(),runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        p=s.HERE/'refine65-protocol.json'
        if p.exists():assert s.read(p)==protocol
        else:s.save(p,protocol)
        state=dict(status='running',pid=os.getpid(),planned=4,completed=[],failures={},started_epoch=time.time())
        s.save(s.HERE/'refine65-state.json',state)
        with ProcessPoolExecutor(max_workers=2) as pool:
            fs={pool.submit(s.run,j):j['id'] for j in JOBS}
            for f in as_completed(fs):
                try:f.result();state['completed'].append(fs[f])
                except Exception as e:state['failures'][fs[f]]=repr(e)
                s.save(s.HERE/'refine65-state.json',state)
        state['status']='complete' if not state['failures'] else 'failed';s.save(s.HERE/'refine65-state.json',state)
        if state['failures']:raise SystemExit(1)
    finally:lock.unlink(missing_ok=True)
if __name__=='__main__':main()
