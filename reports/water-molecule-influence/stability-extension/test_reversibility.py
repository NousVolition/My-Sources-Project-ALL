"""Distinguish physical time reversal, numerical reversibility, and dissipation.
Explicit illustrative models; not reconstructions of cropped phase sketches.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from test_stability import integrate

HERE=Path(__file__).resolve().parent


def duffing(gamma):
    def f(z):
        x,v=z
        return np.array([v,x-x**3-gamma*v])
    return f


def energy(z):
    x,v=z[...,0],z[...,1]
    return v*v/2-x*x/2+x**4/4


def verlet(z0,T,h,gamma=0):
    z=np.array(z0,dtype=float);values=[z.copy()]
    decay=np.exp(-gamma*h/2)
    for _ in range(round(T/h)):
        x,v=z
        v*=decay
        v+=h*(x-x**3)/2
        x+=h*v
        v+=h*(x-x**3)/2
        z=np.array([x,v*decay]);values.append(z.copy())
    return np.arange(len(values))*h,np.array(values)


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    initial=np.array([.5,.2]);R=np.array([1.,-1.]);rows=[];arrays={}
    for gamma in [0.,.2]:
        for method in ['RK4','Verlet']:
            for h in [.05,.025,.0125]:
                evolve=lambda z:integrate(duffing(gamma),z,5.,h) if method=='RK4' else verlet(z,5.,h,gamma)
                t,forward=evolve(initial)
                _,back=evolve(R*forward[-1])
                returned=R*back[-1]
                defect=float(np.linalg.norm(returned-initial))
                if gamma==0 and method=='Verlet':assert defect<1e-11
                if gamma==.2:assert defect>.1
                key=f'g{gamma:g}_{method}_h{h:g}'
                arrays[key+'_time'],arrays[key+'_forward'],arrays[key+'_back']=t,forward,back
                rows.append({'gamma':gamma,'method':method,'h':h,'round_trip_defect':defect,
                             'forward_energy_change':float(energy(forward[-1])-energy(initial)),
                             'maximum_forward_energy_error':float(np.max(abs(energy(forward)-energy(initial))))})
    rk=[r['round_trip_defect'] for r in rows if r['gamma']==0 and r['method']=='RK4']
    assert rk[1]<rk[0]/10 and rk[2]<rk[1]/10
    # Dissipation identity for the mechanical model: dH/dt=-gamma*v^2.
    def augmented(t,z):
        x,v,loss=z
        return [v,x-x**3-.2*v,.2*v*v]
    reference=solve_ivp(augmented,[0,5],[*initial,0.],method='DOP853',rtol=1e-12,atol=1e-14)
    assert reference.success
    energy_balance=float(energy(reference.y[:2,-1])-energy(initial)+reference.y[2,-1])
    assert abs(energy_balance)<1e-10

    # Reversible but not area preserving: R=-I and F(Rz)=-R F(z).
    def nonconservative(z):
        x,y=z
        return np.array([1-x*x,-x*y])
    for z in [[.2,.5],[-.7,1.2],[1.,0.]]:
        z=np.array(z);assert np.allclose(nonconservative(-z),nonconservative(z))
    z0=np.array([.2,.5]);t,forward=integrate(nonconservative,z0,2.,.005)
    _,returned=integrate(nonconservative,-forward[-1],2.,.005)
    analytic=np.stack([np.tanh(t+np.arctanh(z0[0])),z0[1]*np.cosh(np.arctanh(z0[0]))/np.cosh(t+np.arctanh(z0[0]))],axis=-1)
    nc_reference_error=float(np.max(np.linalg.norm(forward-analytic,axis=1)))
    nc_defect=float(np.linalg.norm(-returned[-1]-z0))
    assert nc_reference_error<1e-8 and nc_defect<1e-7
    arrays['nonconservative_time'],arrays['nonconservative_forward']=t,forward
    def exercise_661(z):
        x,y=z
        return np.array([y*(1-x*x),1-y*y])
    def exercise_662(z):
        x,y=z
        return np.array([y,x*np.cos(y)])
    exercises=[]
    for name,f,z0,T in [('6.6.1',exercise_661,np.array([.2,.3]),1.5),
                         ('6.6.2',exercise_662,np.array([.2,.3]),1.5)]:
        R=np.array([1.,-1.])
        for point in [[.2,.3],[-.7,1.2],[1.5,-1.8]]:
            point=np.array(point)
            assert np.allclose(f(R*point),-R*f(point))
        errors=[]
        for h in [.02,.01,.005]:
            t,z=integrate(f,z0,T,h)
            _,back=integrate(f,R*z[-1],T,h)
            errors.append(float(np.linalg.norm(R*back[-1]-z0)))
        ref=solve_ivp(lambda t,z:f(z),[0,T],z0,t_eval=t,method='DOP853',rtol=1e-12,atol=1e-14)
        reference_error=float(np.max(np.linalg.norm(ref.y.T-z,axis=1)))
        assert ref.success and errors[-1]<1e-7 and reference_error<1e-8
        record={'exercise':name,'reversal':'R(x,y)=(x,-y)','h':[.02,.01,.005],
                'round_trip_defects':errors,'RK4_DOP853_max_distance':reference_error}
        if name=='6.6.1':
            record['fixed_points']=[{'point':[1,1],'eigenvalues':[-2,-2],'type':'stable node'},
                                    {'point':[-1,1],'eigenvalues':[2,-2],'type':'saddle'},
                                    {'point':[1,-1],'eigenvalues':[2,2],'type':'unstable node'},
                                    {'point':[-1,-1],'eigenvalues':[-2,2],'type':'saddle'}]
            c=np.arctanh(z0[1])
            exact=np.stack([np.tanh(np.arctanh(z0[0])+np.log(np.cosh(t+c)/np.cosh(c))),np.tanh(t+c)],axis=1)
            record['RK4_analytic_max_distance']=float(np.max(np.linalg.norm(z-exact,axis=1)))
            assert record['RK4_analytic_max_distance']<1e-8
        else:
            record['fixed_points']=[{'point':[0,0],'eigenvalues':[1,-1],'type':'saddle'}]
            record['invariant_lines']='y=(k+1/2)*pi, where x_dot=y and y_dot=0.'
        arrays[name+'_time'],arrays[name+'_forward']=t,z
        exercises.append(record)
    result={'mechanical_model':'x_dot=v; v_dot=x-x^3-gamma*v',
            'mechanical_energy':'H=v^2/2-x^2/2+x^4/4; dH/dt=-gamma*v^2',
            'protocol':'Initial (x,v)=(0.5,0.2); evolve for 5, flip v, evolve forward for 5 with the same gamma, flip v again, compare with the initial state.',
            'rows':rows,'independent_energy_balance_residual':energy_balance,
            'reversible_nonconservative':{'equation':'x_dot=1-x^2, y_dot=-x*y',
              'reversal':'R(x,y)=(-x,-y)','fixed_points':[[-1,0],[1,0]],
              'types':['unstable node','stable node'],'divergence':'-3*x',
              'divergence_at_fixed_points':[3,-3],
              'analytic_solution_max_error':nc_reference_error,'round_trip_defect':nc_defect},
            'interpretation':'Time-reversal symmetry depends on the specified involution R. It does not imply area preservation. A symmetric integrator can reverse accurately without preserving exact energy. Damping breaks velocity-flip reversibility in the mechanical model.',
            'supplied_exercises':exercises,
            'source':'Exercises 6.6.1 and 6.6.2 are transcribed from the supplied page. The mechanical and nonconservative auxiliary models are explicitly chosen benchmarks; other unidentified cropped systems are not reconstructed.',
            'passed':True}
    np.savez_compressed(out/'reversibility.npz',**arrays)
    (out/'reversibility.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run(HERE/'data')
