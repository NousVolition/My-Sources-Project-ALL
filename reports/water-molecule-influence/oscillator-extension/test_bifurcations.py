"""Explicit normal-form controls, not reconstructions of the cropped bifurcation sketch."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from test_oscillators import checked_solve

HERE=Path(__file__).resolve().parent


def run():
    out=HERE/'data';out.mkdir(parents=True,exist_ok=True)
    arrays={};rows=[];pitchfork=[];green=[]
    for mu in [-.2,0.,.2]:
        for r0 in [.1,1.]:
            def f(t,z):
                x,y=z;r2=x*x+y*y
                return [mu*x-y-x*r2,x+mu*y-y*r2]
            t=np.linspace(0,80,8001)
            sol=checked_solve(f,[0,80],[r0,0.],method='DOP853',rtol=1e-11,atol=1e-13,
                              max_step=.1,t_eval=t)
            if mu==0:
                radius=r0/np.sqrt(1+2*r0*r0*t)
            else:
                radius=np.sqrt(mu/(1+(mu/r0**2-1)*np.exp(-2*mu*t)))
            exact=np.stack([radius*np.cos(t),radius*np.sin(t)],axis=1)
            error=float(np.max(np.linalg.norm(sol.y.T-exact,axis=1)))
            assert error<1e-8
            if mu>0:assert abs(np.linalg.norm(sol.y[:,-1])-np.sqrt(mu))<1e-8
            rows.append({'mu':mu,'initial_radius':r0,'final_radius':float(np.linalg.norm(sol.y[:,-1])),
                         'exact_max_distance':error,'origin_trace':2*mu,'origin_determinant':mu*mu+1,
                         'origin_eigenvalues':[[mu,1.],[mu,-1.]],
                         'stable_cycle_radius':float(np.sqrt(mu)) if mu>0 else None})
            arrays[f'hopf_{mu:g}_{r0:g}_t'],arrays[f'hopf_{mu:g}_{r0:g}_z']=t,sol.y.T
        points=[0.] if mu<=0 else [-np.sqrt(mu),0.,np.sqrt(mu)]
        pitchfork.append({'mu':mu,'equilibria':[
            {'point':[float(x),0.],'eigenvalues':[float(mu-3*x*x),-1.]}
            for x in points]})
    # Green's theorem under the explicitly approximate circular-cycle assumption.
    for radius in [1.,2.,3.]:
        numerical=quad(lambda x:2*np.sqrt(max(0.,radius*radius-x*x))*(1-x*x),
                       -radius,radius,epsabs=1e-10,epsrel=1e-12)[0]
        formula=np.pi*radius*radius*(1-radius*radius/4)
        assert abs(numerical-formula)<1e-8
        green.append({'disk_radius':radius,'integral_divergence_divided_by_mu':float(numerical),
                      'analytic_formula':float(formula)})
    result={'hopf':{'model':'x_dot=mu*x-y-x*(x^2+y^2); y_dot=x+mu*y-y*(x^2+y^2)',
                   'polar':'r_dot=mu*r-r^3; theta_dot=1','initial_radii':[.1,1.],
                   'time_end':80.,'solver':'DOP853; rtol=1e-11, atol=1e-13, max_step=0.1',
                   'rows':rows,'interpretation':'At mu=0 the eigenvalues are +/-i, not zero; determinant=1. For mu>0 the origin is unstable and the stable periodic orbit has radius sqrt(mu). At mu=0 cubic damping makes the origin asymptotically stable.'},
            'pitchfork':{'model':'x_dot=mu*x-x^3; y_dot=-y','rows':pitchfork,
                         'interpretation':'At mu=0 the origin has eigenvalues 0,-1. For positive mu two stable equilibria appear and the origin becomes a saddle.'},
            'green_theorem':{'model':'van der Pol x_dot=v, v_dot=mu*(1-x^2)*v-x',
                            'divergence':'mu*(1-x^2)',
                            'circular_region_integral':'mu*pi*r^2*(1-r^2/4)',
                            'nonzero_predicted_radius':2.,'quadrature':green,
                            'limitation':'The circular boundary is a weak-nonlinearity approximation; finite-mu cycles are not exactly circular. Zero net flux is necessary for the actual cycle, not sufficient to prove existence or attraction.'},
            'source':'Explicit illustrative normal forms chosen to test the distinction described in the supplied bifurcation page; the cropped phase portrait has no complete equation.',
            'passed':True}
    (out/'bifurcations.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(out/'bifurcations.npz',**arrays)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':run()
