"""Run only this frozen, bounded comparison; never overwrite results."""
from pathlib import Path
import hashlib,json,sys,time
import model
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parent
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(name,obj):(ROOT/'results'/name).write_text(json.dumps(obj,indent=2),encoding='utf-8')
def main():
    freeze=json.loads((ROOT/'freeze.json').read_text())
    for name,sha in freeze['files'].items():assert digest(ROOT/name)==sha,name
    p=json.loads((ROOT/'protocol.json').read_text())
    (ROOT/'results').mkdir(exist_ok=False)
    start=time.monotonic();cases=model.make_cases(p);write('cases.json',cases)
    assert len(cases)==p['expected']['total']
    t=np.arange(round(p['time']['end']/p['time']['output_dt'])+1)*p['time']['output_dt']
    histories=[];stage_records=[]
    for dt in p['time']['steps']:
        y,record=model.integrate(p,cases,dt)
        histories.append(y);stage_records.append(record)
        print('Integrated',len(cases),'trajectories at dt',dt,flush=True)
    fine=histories[-1];tol=p['verification'];checks=[]
    for a,b in zip(range(len(histories)-1),range(1,len(histories))):
        err=np.max(np.abs(histories[a]-histories[b]),axis=(0,2))
        checks.append(dict(name='step_'+str(p['time']['steps'][a])+'_to_'+str(p['time']['steps'][b]),
                           maximum=float(err.max()),limit=tol['half_step_max_abs'],passed=bool(np.all(err<=tol['half_step_max_abs'])),
                           per_case=err.tolist()))
    independent=[]
    for case in cases:
        if case['pair'] not in (0,6):continue
        if time.monotonic()-start>p['resources']['max_elapsed_seconds']:
            raise TimeoutError('Whole-study wall budget')
        sol=solve_ivp(model.independent_rhs(p,case),(0,p['time']['end']),case['x0']+p['guard_initial'],
                      method='DOP853',t_eval=t,rtol=1e-11,atol=1e-13,max_step=.025)
        err=float(np.max(np.abs(sol.y.T-fine[:,case['id']]))) if sol.success else float('inf')
        independent.append(dict(id=case['id'],success=bool(sol.success),error=err,passed=bool(sol.success and err<=tol['independent_max_abs'])))
    checks.append(dict(name='independent_DOP853',maximum=max(c['error'] for c in independent),
                       limit=tol['independent_max_abs'],passed=all(c['passed'] for c in independent),records=independent))
    assert len(independent)==p['expected']['independent_solver_checks']
    lookup={(c['pair'],c['mode'],c['weight'],c['disturbance']):c['id'] for c in cases}
    swap_weight={'equal':'equal','B_to_A_stronger':'A_to_B_stronger','A_to_B_stronger':'B_to_A_stronger'}
    relabel=[]
    for c in cases:
        if c['pair']>=6:continue
        pulse=c['disturbance']
        if pulse!='none':pulse=('B' if pulse.startswith('A') else 'A')+pulse[1:]
        other=lookup[(c['pair']+6,c['mode'],swap_weight[c['weight']],pulse)]
        relabel.append(float(np.max(np.abs(fine[:,c['id']][:,[1,0,3,2]]-fine[:,other]))))
    checks.append(dict(name='participant_relabel',maximum=max(relabel),limit=tol['relabel_max_abs'],passed=max(relabel)<=tol['relabel_max_abs']))
    shared=[c['id'] for c in cases if c['mode']=='shared']
    equal=[c['id'] for c in cases if c['mode']=='shared' and c['weight']=='equal']
    equal_err=float(np.max(np.abs(fine[:,equal,2]-fine[:,equal,3])))
    checks.append(dict(name='equal_realized_gains_shared_equal_weights',maximum=equal_err,limit=tol['equal_transmission_tolerance'],passed=equal_err<=tol['equal_transmission_tolerance']))
    shared_err=float(np.max(np.abs(fine[:,shared,2]-fine[:,shared,3])))
    checks.append(dict(name='shared_guard_remains_shared',maximum=shared_err,limit=tol['equal_transmission_tolerance'],passed=shared_err<=tol['equal_transmission_tolerance']))
    analytical=[]
    for c in cases:
        if c['mode']!='frozen' or c['disturbance']!='none':continue
        if time.monotonic()-start>p['resources']['max_elapsed_seconds']:
            raise TimeoutError('Whole-study wall budget')
        ca,cb=np.array(c['k'])*p['guard_initial']
        lam=p['lambda'];M=np.array([[-lam-ca,ca],[cb,-lam-cb]])
        eq=np.linalg.solve(-M,lam*np.array(c['b']))
        exact=np.array([eq+expm(M*tt)@(np.array(c['x0'])-eq) for tt in t])
        analytical.append(float(np.max(np.abs(exact-fine[:,c['id'],:2]))))
    checks.append(dict(name='frozen_guard_exact_matrix_exponential',maximum=max(analytical),limit=tol['analytic_max_abs'],passed=max(analytical)<=tol['analytic_max_abs']))
    pulse_t=np.linspace(p['pulse']['onset'],p['pulse']['onset']+p['pulse']['duration'],10001)
    pulse=p['pulse']['amplitude']*np.sin(np.pi*(pulse_t-p['pulse']['onset'])/p['pulse']['duration'])**2
    checks.append(dict(name='pulse_area',value=float(np.trapezoid(pulse,pulse_t)),passed=abs(float(np.trapezoid(pulse,pulse_t))-p['pulse']['signed_area_magnitude'])<1e-12))
    np.savez_compressed(ROOT/'results/trajectories.npz',time=t,state=fine,
                        coarse=histories[0],half=histories[1])
    write('verification.json',dict(passed=all(c['passed'] for c in checks),checks=checks,stages=stage_records,elapsed_seconds=time.monotonic()-start,
                                   protocol_sha256=digest(ROOT/'protocol.json'),frozen_files={n:digest(ROOT/n) for n in freeze['files']}))
    print('Verification:',all(c['passed'] for c in checks),'wall seconds',time.monotonic()-start,flush=True)
    return 0 if all(c['passed'] for c in checks) else 2

if __name__=='__main__':
    try:sys.exit(main())
    except Exception as exc:
        if (ROOT/'results').exists():write('failure.json',{'type':type(exc).__name__,'error':str(exc)})
        raise

