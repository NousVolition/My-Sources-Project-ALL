"""Retain the initial failed cross-checks; add tighter, independently checked runs."""
from pathlib import Path
import json
import numpy as np
from scipy.special import ellipj, ellipkinc
from stress_solver import model, resolved, scaled_error, energy_budget

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'

def conservative_exact(t,q0=2.,v0=3.):
    # q'' + q + q^3 = 0.  Jacobi cn has this cubic restoring law.
    e=.5*v0*v0+.5*q0*q0+.25*q0**4
    amplitude=np.sqrt(-1+np.sqrt(1+4*e))
    omega=np.sqrt(1+amplitude**2)
    modulus=amplitude**2/(2*omega**2)
    phase=-ellipkinc(np.arccos(q0/amplitude),modulus)
    sn,cn,dn,_=ellipj(omega*np.asarray(t)+phase,modulus)
    y=np.zeros((len(t),5));y[:,0]=amplitude*cn;y[:,1]=-amplitude*omega*sn*dn
    return y

def main():
    findings=[]
    for name in ('long-conservative','long-weak-damping'):
        row=json.loads((DATA/(name+'.json')).read_text())
        old=np.load(DATA/(name+'.npz'));t=old['reference_time']
        p=model.Parameters(**row['parameters']);ys={};infos={}
        for method in ('DOP853','Radau'):
            evaluate,info=resolved(p,row['q0'],row['v0'],row['end'],method=method,rtol=2e-12)
            ys[method]=evaluate(t);infos[method]=info
            np.savez_compressed(DATA/(name+'-tight-'+method+'.npz'),time=t,state=ys[method])
            print(name,method,'completed',flush=True)
        error=scaled_error(ys['Radau'],ys['DOP853'])
        item=dict(name=name,threshold=2e-6,scaled_discrepancy=error,passed=error<2e-6,
            solvers=infos,original_reference_change=scaled_error(old['reference_state'],ys['DOP853']),
            energy={k:energy_budget(y,p) for k,y in ys.items()})
        # Compare the original RK4 trajectories at every saved step with the tighter DOP853.
        item['rk4_errors']={}
        for tag in ('coarse','fine'):
            tt=old[tag+'_time'];ref=evaluate(tt)  # final evaluate is the tight Radau solution
            item['rk4_errors'][tag]=scaled_error(old[tag+'_state'],ref)
        if name=='long-conservative':
            exact=conservative_exact(t)
            np.savez_compressed(DATA/(name+'-exact.npz'),time=t,state=exact)
            item['exact_errors']={k:scaled_error(y,exact) for k,y in ys.items()}
            item['original_exact_error']=scaled_error(old['reference_state'],exact)
            item['exact_initial_error']=float(max(abs(exact[0]-[2,3,0,0,0])))
            item['passed']=item['passed'] and max(item['exact_errors'].values())<2e-6
        findings.append(item)
    (DATA/'long-refinement.json').write_text(json.dumps(findings,indent=2)+'\n')
    print(json.dumps(findings,indent=2))
    if not all(x['passed'] for x in findings):raise SystemExit(1)

if __name__=='__main__':main()
