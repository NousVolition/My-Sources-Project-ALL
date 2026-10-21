from pathlib import Path
import os,json,hashlib,time,traceback
from datetime import datetime,timezone
import numpy as np
from scipy.integrate import solve_ivp
from model import Model,cases,run
R=Path(__file__).resolve().parent
def serial(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x))
def write(path,x):
    path=R/path;temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(x,indent=2,default=serial),encoding='utf-8');temp.replace(path)
def verify():
    f=json.loads((R/'freeze.json').read_text())
    for name,sha in f['files'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==sha,name
def main():
    verify();p=json.loads((R/'protocol.json').read_text());cs=cases(p)
    assert len(cs)==45
    (R/'results').mkdir(exist_ok=False)
    with (R/'controller.lock').open('x',encoding='utf-8') as f:json.dump(dict(pid=os.getpid(),created_utc=datetime.now(timezone.utc).isoformat()),f)
    write('handoff.json',dict(status='running',pid=os.getpid(),protocol='protocol.json',freeze='freeze.json',summary='results/summary.json',stdout='run.log',stderr='run.err',no_automatic_retry=True))
    write('results/cases.json',cs);start=time.monotonic();completed=[]
    for batch in p['batches']:
        n=batch['n'];dt=batch['dt'];label=f'n{n}-dt{dt:g}'
        def progress(t,bounds):write('results/progress.json',dict(status='running',batch=label,t=t,pid=os.getpid(),completed=completed,bounds=bounds))
        progress(0,{})
        d,fields,times=run(p,cs,n,dt,progress)
        np.savez_compressed(R/'results'/f'{label}.npz',fields=fields,time=times)
        write(f'results/{label}.json',d)
        completed.append(label);print(label,'complete',d['elapsed_seconds'],'seconds',flush=True)
        if time.monotonic()-start>p['resources']['max_total_seconds']:raise TimeoutError('Total study budget')
    # Relabel the same directional cases; it is a control, not another sample.
    original=[c for c in cs if c['mode']=='directional']
    swapped=[dict(c,k=c['k'][::-1],pulse={'none':'none','A':'B','B':'A'}[c['pulse']]) for c in original]
    write('results/progress.json',dict(status='running',batch='participant_relabel',pid=os.getpid(),completed=completed))
    d,fields,times=run(p,swapped,32,.01,swapped=True)
    np.savez_compressed(R/'results/relabel.npz',fields=fields,time=times)
    write('results/relabel.json',d)
    # Independent evolution: different momentum form, scalar transport form,
    # time integrator and adaptive step selection.
    independent=[];indfields=[];times=np.array(p['time']['field_times'])
    reference=np.load(R/'results/n24-dt0.01.npz')['fields']
    for c in [c for c in original if c['pulse']=='B']:
        write('results/progress.json',dict(status='running',batch='independent',case=c['id'],pid=os.getpid(),completed=completed))
        m=Model(24,p,[c]);h=m.start();shape=h.shape
        sol=solve_ivp(lambda t,z:m.independent_rhs(t,z.reshape(shape)).ravel(),(0,p['time']['end']),h.ravel(),
            method='DOP853',t_eval=times,rtol=p['independent']['rtol'],atol=p['independent']['atol'],max_step=p['independent']['max_step'])
        if sol.success:
            ff=np.array([m.real(z.reshape(shape))[0] for z in sol.y.T]);err=float(np.max(np.abs(ff-reference[:,c['id']])))
        else:ff=np.full_like(reference[:,c['id']],np.nan);err=float('inf')
        independent.append(dict(case=c['id'],success=bool(sol.success),message=sol.message,max_field_error=err,passed=bool(sol.success and err<p['independent']['limit']),evaluations=sol.nfev))
        indfields.append(ff)
    np.savez_compressed(R/'results/independent.npz',fields=np.array(indfields),time=times)
    write('results/independent.json',independent);verify()
    summary=dict(status='completed',batch_labels=completed,conditions=len(cs),elapsed_seconds=time.monotonic()-start,
                 independent_checks=independent,interpretation='Await reporting and all numerical/scientific screens; completion alone is not verification.')
    write('results/summary.json',summary);write('results/progress.json',dict(status='completed',pid=os.getpid(),completed=completed))
    write('handoff.json',dict(status='integration_completed',pid=os.getpid(),protocol='protocol.json',freeze='freeze.json',summary='results/summary.json',stdout='run.log',stderr='run.err',next='Analyze saved results; retain all failed screens. Do not rerun integrations.'))
    print('COMPLETE',summary['elapsed_seconds'],flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        error=traceback.format_exc()
        if (R/'results').is_dir():write('results/failure.json',dict(error=error,pid=os.getpid()))
        print(error,flush=True);raise
