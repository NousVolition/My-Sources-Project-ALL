"""RK4 stability test: exact linear calibration, nonlinear local limits,
and rabbit/sheep basin prediction. All parameters and starting states are fixed.
"""
from pathlib import Path
import argparse,json,platform
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
SYSTEMS=[
 ('5.2.3',[[0,1],[-2,-3]],[-1,-2],'asymptotically stable node',[(1,-1),(1,-2)]),
 ('5.2.4',[[5,10],[-1,-1]],[2+1j,2-1j],'unstable spiral',[]),
 ('5.2.5',[[3,-4],[1,-1]],[1,1],'unstable defective node',[(2,1)]),
 ('5.2.6',[[-3,2],[1,-2]],[-1,-4],'asymptotically stable node',[(1,1),(2,-1)]),
 ('5.2.7',[[5,2],[-17,-5]],[3j,-3j],'center: Lyapunov stable, not attracting',[]),
 ('5.2.8',[[-3,4],[-2,3]],[1,-1],'saddle',[(1,1),(2,1)]),
 ('5.2.9',[[4,-3],[8,-6]],[0,-2],'stable line of equilibria; individual points not attracting',[(3,4),(1,2)]),
 ('5.2.10',[[1,0],[1,-2]],[1,-2],'saddle',[(3,1),(0,1)])]


def rk4_step(f,x,h):
    k1=h*f(x);k2=h*f(x+k1/2);k3=h*f(x+k2/2);k4=h*f(x+k3)
    return x+(k1+2*k2+2*k3+k4)/6


def integrate(f,x0,T,h,save=True):
    steps=int(round(T/h))
    assert steps>0 and abs(steps*h-T)<1e-10
    x=np.array(x0,dtype=float);values=[x.copy()] if save else None
    for _ in range(steps):
        x=rk4_step(f,x,h)
        if save:values.append(x.copy())
    assert np.isfinite(x).all()
    if save:return np.arange(steps+1)*h,np.array(values)
    return x


def competition(z):
    x,y=z[...,0],z[...,1]
    return np.stack([x*(3-x-2*y),y*(2-x-y)],axis=-1)


def competition_jacobian(x,y):
    return np.array([[3-2*x-2*y,-2*x],[-y,2-x-2*y]])


def winners(z):
    rabbit=np.linalg.norm(z-[3,0],axis=-1)<1e-5
    sheep=np.linalg.norm(z-[0,2],axis=-1)<1e-5
    return np.where(rabbit,1,np.where(sheep,2,0))


def stable_manifold(delta,h=.002):
    v=np.array([np.sqrt(2),1]);v/=np.linalg.norm(v)
    branches=[]
    for sign in [-1,1]:
        z=np.array([1.,1.])+sign*delta*v;values=[z.copy()]
        for _ in range(25000):
            z=rk4_step(competition,z,-h)
            values.append(z.copy())
            if np.max(z)>4.05 or np.linalg.norm(z)<1e-5:break
        else:raise AssertionError('Manifold boundary not reached')
        branch=np.array(values)
        assert np.isfinite(branch).all() and np.min(branch)>-1e-12
        branches.append(branch)
    return branches


def cubic(z):
    x,y=z[...,0],z[...,1]
    return np.stack([-x+x**3,-2*y],axis=-1)


def cubic_exact(z,T):
    x,y=z
    denominator=1-x*x+x*x*np.exp(-2*T)
    assert denominator>0
    return np.array([x*np.exp(-T)/np.sqrt(denominator),y*np.exp(-2*T)])


def radial(a):
    def f(z):
        r2=np.sum(z*z,axis=-1)
        return np.stack([-z[...,1]+a*z[...,0]*r2,z[...,0]+a*z[...,1]*r2],axis=-1)
    return f


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    arrays={};linear=[]
    initials=np.array([[1.,.3],[-.7,1.]])
    for name,raw,eig,classification,vectors in SYSTEMS:
        A=np.array(raw,dtype=float)
        trace,det=float(np.trace(A)),float(np.linalg.det(A))
        assert np.allclose(np.poly(A),np.poly(np.array(eig,dtype=complex)).real)
        for vector in vectors:
            v=np.array(vector,dtype=float);w=A@v
            assert abs(np.linalg.det(np.stack([v,w])))<1e-10
        exact=np.einsum('ij,nj->ni',expm(A),initials)
        errors=[]
        for h in [.05,.025,.0125]:
            final=integrate(lambda z:z@A.T,initials,1.,h,save=False)
            errors.append(float(np.linalg.norm(final-exact)/np.linalg.norm(exact)))
        orders=np.log2(np.array(errors[:-1])/errors[1:])
        assert errors[-1]<1e-5 and np.all((orders>3.5)&(orders<4.5))
        linear.append({'exercise':name,'A':raw,'trace':trace,'determinant':det,
                       'eigenvalues':[str(z) for z in eig],'classification':classification,
                       'real_eigenvectors':vectors,'relative_errors':errors,'observed_orders':orders.tolist()})
    # Defective eigenvalue and the non-isolated equilibrium case.
    A=np.array(SYSTEMS[2][1]);v=np.array([2,1]);g=np.array([1,0])
    assert np.array_equal((A-np.eye(2))@g,v)
    A9=np.array(SYSTEMS[6][1],dtype=float);P=np.eye(2)+A9/2
    assert np.allclose(P@P,P) and np.allclose(A9@P,0)
    A7=np.array(SYSTEMS[4][1],dtype=float);Q=np.array([[17,5],[5,2.]])
    assert np.allclose(A7.T@Q+Q@A7,0) and np.linalg.eigvalsh(Q).min()>0
    energy_controls=[]
    for h in [.1,.05,.025]:
        t,z=integrate(lambda z:z@A7.T,[.2,0.],200,h)
        energy=np.einsum('ni,ij,nj->n',z,Q,z)
        exact_energy=np.array([.2,0.])@Q@np.array([.2,0.])
        energy_controls.append({'h':h,'final_relative_invariant_error':float(energy[-1]/exact_energy-1)})
        arrays[f'center_h{h:g}_time'],arrays[f'center_h{h:g}_energy']=t,energy/exact_energy

    # Supplied Jacobian + equilibria reconstruct f=(-x+x^3,-2y).
    local=[]
    for point,T,radii in [([0.,0.],1.,[.2,.1,.05,.025]),([1.,0.],.25,[.05,.025,.0125,.00625])]:
        J=np.diag([-1+3*point[0]**2,-2.])
        errors=[]
        for rho in radii:
            initial=np.array(point)+[rho,rho]
            exact=cubic_exact(initial,T)
            predicted=np.array(point)+expm(J*T)@np.array([rho,rho])
            numerical=integrate(cubic,initial,T,.001,False)
            error=float(np.linalg.norm(exact-predicted));errors.append(error)
            local.append({'fixed_point':point,'T':T,'rho':rho,'exact_linearization_error':error,
                          'RK4_error':float(np.linalg.norm(numerical-exact))})
        expected=3 if point[0]==0 else 2
        assert abs(np.polyfit(np.log(radii),np.log(errors),1)[0]-expected)<.15
    # Chosen nonlinear completion of a center: all a share J=[[0,-1],[1,0]].
    radial_rows=[]
    for a in [-1.,0.,1.]:
        t,z=integrate(radial(a),[.1,0.],40.,.025)
        exact_r=.1/np.sqrt(1-2*a*.1**2*t)
        error=float(np.max(abs(np.linalg.norm(z,axis=1)-exact_r)))
        assert error<1e-5
        arrays[f'radial_a{a:g}_time'],arrays[f'radial_a{a:g}_radius']=t,np.linalg.norm(z,axis=1)
        radial_rows.append({'a':a,'initial_radius':.1,'time_end':40.,'max_radius_error':error,
                            'final_radius':float(np.linalg.norm(z[-1]))})
    escapes=[]
    for r0 in [.2,.1,.05]:
        exact=(r0**-2-.5**-2)/2
        for h in [.05,.025]:
            z=np.array([r0,0.]);previous=r0;elapsed=0.
            for _ in range(int(np.ceil(1.2*exact/h))+1):
                nxt=rk4_step(radial(1),z,h);current=np.linalg.norm(nxt)
                if current>=.5:
                    estimate=elapsed+h*(.5-previous)/(current-previous)
                    break
                z=nxt;previous=current;elapsed+=h
            else:raise AssertionError('Expected escape not detected')
            assert abs(estimate-exact)<2*h
            escapes.append({'initial_radius':r0,'epsilon':.5,'h':h,'exact_exit_time':exact,
                            'interpolated_exit_time':float(estimate),'crossing_step':[elapsed,elapsed+h]})

    # Full population model: two attractors separated by the saddle's stable manifold.
    equilibrium_rows=[]
    for point in [[0,0],[0,2],[3,0],[1,1]]:
        J=competition_jacobian(*point)
        assert np.allclose(competition(np.array(point,dtype=float)),0)
        numerical_J=np.column_stack([(competition(np.array(point)+1e-6*np.eye(2)[j])-competition(np.array(point)-1e-6*np.eye(2)[j]))/2e-6 for j in range(2)])
        assert np.max(abs(J-numerical_J))<1e-8
        eigenvalues=np.linalg.eigvals(J)
        assert np.max(abs(eigenvalues.imag))<1e-12
        equilibrium_rows.append({'point':point,'J':J.tolist(),'eigenvalues':[float(z.real) for z in eigenvalues]})
    coordinates=np.linspace(.05,3.5,47)
    X,Y=np.meshgrid(coordinates,coordinates);starts=np.stack([X.ravel(),Y.ravel()],axis=-1)
    final_coarse=integrate(competition,starts,80.,.02,False)
    final_fine=integrate(competition,starts,80.,.01,False)
    coarse,fine=winners(final_coarse),winners(final_fine)
    assert np.all(coarse==fine) and np.all(fine>0)
    tangent=np.where(starts[:,1]-1>(starts[:,0]-1)/np.sqrt(2),2,1)
    mismatch=tangent!=fine
    distance=np.linalg.norm(starts-[1,1],axis=1)
    agreement=[]
    for radius in [.25,.5,1.,4.]:
        mask=distance<radius
        agreement.append({'radius_about_saddle':radius,'samples':int(mask.sum()),
                           'tangent_prediction_errors':int(np.count_nonzero(mismatch&mask))})
    manifolds=stable_manifold(1e-5);refined=stable_manifold(1e-6)
    def combine(branches):
        c=np.concatenate([branches[0][::-1],[[1,1]],branches[1]])
        assert np.all(np.diff(c[:,0])>=0)
        return c
    curve,curve2=combine(manifolds),combine(refined)
    probes=np.array([.25,.5,.75,1.5,2.,3.])
    manifold_error=float(np.max(abs(np.interp(probes,curve[:,0],curve[:,1])-np.interp(probes,curve2[:,0],curve2[:,1]))))
    assert manifold_error<1e-4
    ysep=float(np.interp(.5,curve2[:,0],curve2[:,1]))
    near=np.array([[.5,ysep-.001],[.5,ysep+.001]])
    tt,trajectories=integrate(competition,near,80.,.01)
    assert np.array_equal(winners(trajectories[-1]),[1,2])
    checks=[]
    for i,start in enumerate(near):
        independent=solve_ivp(lambda t,z:competition(z),[0,80],start,t_eval=tt,method='DOP853',rtol=1e-11,atol=1e-13)
        assert independent.success
        error=float(np.max(np.linalg.norm(independent.y.T-trajectories[:,i,:],axis=1)))
        assert error<2e-6
        checks.append({'initial':start.tolist(),'RK4_DOP853_max_distance':error,'final':trajectories[-1,i].tolist()})
    arrays.update(grid_starts=starts,grid_final=final_fine,grid_labels=fine,
                  grid_tangent_labels=tangent,manifold=curve2,
                  near_boundary_time=tt,near_boundary_trajectories=trajectories)
    results={'linear':linear,'linear_steps':[.05,.025,.0125],
             'center_invariant':{'Q':Q.tolist(),'euclidean_stability_bound':float(np.sqrt(np.linalg.eigvalsh(Q).max()/np.linalg.eigvalsh(Q).min())),
                                 'RK4_artificial_drift':energy_controls},
             'line_of_equilibria':{'equation':'y=4x/3','limit_projection':P.tolist()},
             'cubic_local_tests':local,'radial':radial_rows,'escape_tests':escapes,
             'competition':{'equilibria':equilibrium_rows,'grid_count':len(starts),
                            'rabbit_winners':int(np.count_nonzero(fine==1)),'sheep_winners':int(np.count_nonzero(fine==2)),
                            'unresolved':int(np.count_nonzero(fine==0)),'step_halving_outcome_changes':int(np.count_nonzero(coarse!=fine)),
                            'tangent_comparisons':agreement,'manifold_start_refinement_max_y_difference':manifold_error,
                            'separatrix_y_at_x_05':ysep,'near_boundary_checks':checks},'passed':True}
    protocol={'purpose':'Test what local linearization can predict, and identify false numerical evidence of attraction or stability.',
              'integration':'Classical vector RK4 with k_i containing h, as in the supplied page; independent matrix exponential and DOP853 checks.',
              'linear_source':'All eight visible exercises 5.2.3 through 5.2.10, transcribed into the matrices in results.json.',
              'cubic_source':'f=(-x+x^3,-2y), reconstructed from the supplied Jacobian and stated equilibria.',
              'radial_source':'Explicit illustrative completion, not a claimed transcription of the cropped radial figure: x_dot=-y+a*x*(x^2+y^2), y_dot=x+a*y*(x^2+y^2).',
              'radial_exact':'r_dot=a*r^3, theta_dot=1; r(t)=r0/sqrt(1-2*a*r0^2*t). All three Jacobians at the origin have eigenvalues +/-i.',
              'competition_source':'Supplied model x_dot=x*(3-x-2*y), y_dot=y*(2-x-y), x,y>=0.',
              'competition_grid':'47x47 inclusive grid on [0.05,3.5]^2; T=80; RK4 h=0.02 and 0.01. Endpoint tolerance 1e-5 to a stable equilibrium, otherwise unresolved.',
              'separatrix':'Backward RK4 h=0.002 from saddle +/-delta*(sqrt(2),1)/sqrt(3), delta=1e-5 and 1e-6; stop at box boundary 4.05 or distance to origin below 1e-5.',
              'tangent_predictor':'Predict sheep above y=1+(x-1)/sqrt(2), rabbits below. This is a deliberately local approximation tested against the global nonlinear basin.',
              'units':'Illustrative dimensionless state/time. No biological rates fitted. No new water MD.',
              'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
              'randomness':'None; every initial condition and tolerance specified.'}
    np.savez_compressed(out/'trajectories.npz',**arrays)
    for filename,obj in [('results.json',results),('protocol.json',protocol)]:
        (out/filename).write_text(json.dumps(obj,indent=2),encoding='utf-8')
    print(json.dumps({'competition':results['competition'],'radial':radial_rows,'passed':True},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'data')
    args=p.parse_args();run(args.out)
