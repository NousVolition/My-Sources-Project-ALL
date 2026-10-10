"""Independent-start, independent-solver checks, plus exact-cycle controls."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from model import rhs, integrate, voltage, synchronized_voltage
ROOT=Path(__file__).resolve().parent


def main():
    cfg=json.loads((ROOT/'protocol.json').read_text()); c=cfg['reference_validation']
    out=ROOT/'results';out.mkdir(exist_ok=True)
    rng=np.random.default_rng(c['seed']); ps=[];bs=[];aa=[]
    for a in cfg['alphas']:
        for b in c['biases']:
            for _ in range(c['n_starts_per_alpha_bias']):
                ps.append(rng.uniform(-np.pi,np.pi,2));bs.append(b);aa.append(a)
    ps=np.array(ps);bs=np.array(bs);aa=np.array(aa)
    _,p,d=integrate(ps,bs,aa,400,.025)
    fixed=voltage(p[200],p[400],200)
    reference=[];rows=[]
    for j,(p0,b,a) in enumerate(zip(ps,bs,aa)):
        sol=solve_ivp(lambda t,y:rhs(y,b,a),(0,400),p0,method='DOP853',
                      rtol=c['rtol'],atol=c['atol'],max_step=c['max_step'],t_eval=[200,400])
        if not sol.success:raise RuntimeError(sol.message)
        v=voltage(sol.y[:,0],sol.y[:,1],200);reference.append(sol.y.T)
        rows.append({'alpha':float(a),'bias':float(b),'initial':p0.tolist(),
                     'fixed_voltage':fixed[j].tolist(),'reference_voltage':v.tolist(),
                     'total_voltage_error':float(abs((fixed[j]-v).sum())),
                     'max_individual_voltage_error':float(max(abs(fixed[j]-v)))})
    cycles=[]
    for a in cfg['alphas']:
        for b in [1.005,1.2,2.]:
            period=2*np.pi*(1+2*a)/np.sqrt(b*b-1)
            n=int(np.ceil(period/.025))
            _,q,_=integrate(np.zeros(2),b,a,period,period/n,stride=period)
            measured=voltage(q[0],q[-1],period).sum();exact=float(synchronized_voltage(b,a))
            cycles.append({'alpha':a,'bias':b,'period':float(period),'measured':float(measured),
                           'exact':exact,'relative_error':float(abs(measured-exact)/exact)})
    result={'independent_solver_runs':rows,'analytic_cycles':cycles,
            'max_reference_voltage_error':max(max(r['total_voltage_error'],r['max_individual_voltage_error']) for r in rows),
            'max_analytic_cycle_relative_error':max(r['relative_error'] for r in cycles)}
    (out/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/'data/validation.npz',initial=ps,bias=bs,alpha=aa,time=np.arange(401),
                        fixed_phases=p,fixed_dissipated=d,reference_phases=np.array(reference))
    print(json.dumps({k:v for k,v in result.items() if not isinstance(v,list)},indent=2))


if __name__=='__main__':main()
