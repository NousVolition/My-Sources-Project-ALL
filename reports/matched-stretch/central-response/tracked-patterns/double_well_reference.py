"""The user's displayed two-dimensional equation, as a separate RK4 benchmark."""
from pathlib import Path
import os,json
import numpy as np
ROOT=Path(__file__).resolve().parent/'double-well-reference';ROOT.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def force(z):return np.array([z[1],z[0]-z[0]**3])
def energy(z):return .5*z[...,1]**2-.5*z[...,0]**2+.25*z[...,0]**4
def integrate(z0,dt,end=20):
    n=round(end/dt);z=np.array(z0,dtype=float);out=[z.copy()]
    for _ in range(n):
        a=force(z);b=force(z+dt*a/2);c=force(z+dt*b/2);d=force(z+dt*c);z=z+dt*(a+2*b+2*c+d)/6;out.append(z.copy())
    return np.linspace(0,end,n+1),np.array(out)
def main():
    starts=[(-1.2,0),(1.2,0),(1,.6),(0,.5)];records=[]
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    x=np.linspace(-1.9,1.9,25);y=np.linspace(-1.8,1.8,23);xx,yy=np.meshgrid(x,y);vx=yy;vy=xx-xx**3;norm=np.hypot(vx,vy);mask=norm>1e-12
    axs[0].quiver(xx[mask],yy[mask],vx[mask]/norm[mask],vy[mask]/norm[mask],color='.8',pivot='mid',scale=45,width=.002)
    for z0 in starts:
        t,z=integrate(z0,.01);tf,zf=integrate(z0,.005);ef=energy(zf);err=float(np.max(np.linalg.norm(z-zf[::2],axis=1)))
        rec={'start':z0,'initial_energy':float(energy(np.array(z0))),'max_absolute_energy_drift':float(abs(ef-ef[0]).max()),'step_halving_max_state_difference':err,'x_range':[float(zf[:,0].min()),float(zf[:,0].max())]};records.append(rec)
        _,reverse=integrate([zf[-1,0],-zf[-1,1]],.005)
        rec['velocity_reversal_return_error']=float(np.linalg.norm(reverse[-1]-[z0[0],-z0[1]]));assert rec['velocity_reversal_return_error']<1e-7
        label=f'H = {ef[0]:.4f}';axs[0].plot(zf[:,0],zf[:,1],lw=1.4,label=label);axs[1].plot(tf,ef-ef[0],label=label)
        assert rec['max_absolute_energy_drift']<1e-7
        if rec['initial_energy']<0:assert np.all(zf[:,0]*zf[0,0]>0)
        else:assert rec['x_range'][0]<-1 and rec['x_range'][1]>1
    # The homoclinic orbit has an exact time parametrization; finite plots truncate its tails.
    tt=np.linspace(-6,6,2401);xh=np.sqrt(2)/np.cosh(tt);yh=-np.sqrt(2)*np.sinh(tt)/np.cosh(tt)**2
    for sign in (-1,1):axs[0].plot(sign*xh,sign*yh,'k--',lw=1,label='H = 0 separatrix' if sign==1 else None)
    ht,hz=integrate([np.sqrt(2),0],.005,6);exact=np.stack([np.sqrt(2)/np.cosh(ht),-np.sqrt(2)*np.sinh(ht)/np.cosh(ht)**2],axis=1);herr=float(np.max(np.linalg.norm(hz-exact,axis=1)))
    assert herr<1e-6
    axs[0].scatter([0],[0],s=45,marker='x',c='k');axs[0].scatter([-1,1],[0,0],s=32,facecolors='none',edgecolors='k');axs[0].set(xlabel='Position x',ylabel='Velocity y',title='Closed paths and the saddle boundary',xlim=(-1.9,1.9),ylim=(-1.8,1.8));axs[0].legend(fontsize=8,loc='lower left')
    axs[1].set(xlabel='Time',ylabel='Energy change from initial value',title='RK4 energy check, step 0.005');axs[1].legend(fontsize=8);axs[1].grid(alpha=.2)
    fig.suptitle('Double-well reference: dx/dt = y; dy/dt = x - x cubed')
    fig.savefig(ROOT/'phase-plane.png',dpi=155);plt.close(fig)
    result={'equations':{'dx_dt':'y','dy_dt':'x-x^3'},'role':'Separate textbook reference, not the fluid equation','fluid_source_changed':False,'runs':records,'homoclinic_exact_path_max_error_through_6':herr,'fixed_points':[{'point':[0,0],'eigenvalues':[-1,1],'type':'saddle'},{'point':[-1,0],'eigenvalues':['+i sqrt(2)','-i sqrt(2)'],'type':'center by conserved energy'},{'point':[1,0],'eigenvalues':['+i sqrt(2)','-i sqrt(2)'],'type':'center by conserved energy'}]}
    (ROOT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    text=['# Double-well reference test','','**RK4 reproduces the expected separate loops below the energy barrier and a loop around both wells above it. Energy drift and time-step differences are small over the tested interval.**','','This is the exact two-dimensional equation supplied in the screenshot: `dx/dt = y`, `dy/dt = x-x^3`. It is a separate mathematical reference. It does not replace the unforced Navier–Stokes model or add a force to the fluid.','','![Phase plane and energy check](phase-plane.png)','','## Results','','| Start (x,y) | Initial energy | Largest energy drift | Maximum change after halving the RK4 step |','| --- | ---: | ---: | ---: |']
    for r in records:text.append(f"| {tuple(r['start'])} | {r['initial_energy']:.4f} | {r['max_absolute_energy_drift']:.3g} | {r['step_halving_max_state_difference']:.3g} |")
    text+=['','Each regular trajectory runs to time 20 with steps 0.01 and 0.005. Energy is `H = y^2/2 - x^2/2 + x^4/4`; H is used here to avoid confusion with the fluid mirror error E.','','The origin is a saddle, with stable tangent y=-x and unstable tangent y=x. The nonlinear manifolds lie on `H=0`, or `y^2=x^2-x^4/2`. Their two homoclinic loops approach the origin as time tends to positive or negative infinity. The dashed lines use the exact paths `x(t)=+/-sqrt(2) sech(t)`, `y(t)=d x/dt`; the plotted tails are truncated. The maximum numerical difference from that exact path through time 6 is '+f'{herr:.3g}'+'. A finite numerical run cannot complete a homoclinic orbit.','','The points (-1,0) and (1,0) are centers. Conservation of H establishes the nearby closed contours; purely imaginary eigenvalues alone would not establish nonlinear centers. For -1/4 < H < 0, trajectories stay around one well. H > 0 gives an outer orbit encircling both wells. The centers are not attracting resting states.','','## How this helps the fluid study','','This reference checks the distinction between position and velocity, closed motion and settling, energy barriers and attraction, and linear tangent directions versus curved nonlinear manifolds. The fluid measurements can be examined for these patterns, but resemblance of a projection does not establish a double-well potential, a conserved H, or a pair of fluid attractors. No such mapping is fitted here.','','[Computed results](results.json) · [Tracked-fluid analysis](../README.md)','']
    text+=['## Velocity-reversal check','','After each time-20 endpoint, the velocity was negated and the same equation was evolved for another 20 time units at step 0.005. The expected return is the original position with the opposite original velocity. The largest return error was '+f"{max(r['velocity_reversal_return_error'] for r in records):.3g}"+'. This checks the time-reversal symmetry of this undamped reference. It does not establish time reversal for viscous fluid flow, or equate time reversal with a state-dependent change of stretching sign.','']
    (ROOT/'README.md').write_text('\n'.join(text),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
