"""Supplement requested via textbook excerpts: triangular locking and linear dynamics.

All numerical controls here are separate finite-dimensional systems, not NS.
The circuit's RSJ equations are conditional on stated component assumptions.
"""
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.linalg import expm
from core import ROOT,save_json


def triangle(phi):
    q=(np.asarray(phi)+np.pi/2)%(2*np.pi)-np.pi/2
    return np.where(q<=np.pi/2,q,np.pi-q)


def tri_evolve(mu,phi,times,normalized=False):
    scale=2/np.pi if normalized else 1.
    return solve_ivp(lambda t,x:mu-scale*triangle(x),(times[0],times[-1]),np.atleast_1d(phi).astype(float),t_eval=times,method='DOP853',rtol=1e-10,atol=1e-12,max_step=.1).y.T


def classify(A):
    tau=float(np.trace(A));delta=float(np.linalg.det(A));disc=tau*tau-4*delta
    if abs(delta)<1e-12:return 'zero eigenvalue / degenerate'
    if delta<0:return 'saddle'
    if abs(tau)<1e-12 and disc<0:return 'center (linear system)'
    sign='stable' if tau<0 else 'unstable'
    if abs(disc)<1e-12:return sign+' repeated-eigenvalue case'
    return sign+(' node' if disc>0 else ' spiral')


def circuit_rhs(phi,bias=1.2,alpha=.5):
    """Identical overdamped RSJ pair, alpha=r/R, time scaled by 2e Ic r/hbar."""
    M=np.array([[1+alpha,alpha],[alpha,1+alpha]])
    return np.linalg.solve(M,bias-np.sin(phi))


def main():
    out=ROOT/'results';out.mkdir(exist_ok=True);ts=np.linspace(0,12,481)
    tri=[];tracks={}
    for mu in [0.,.5,1.4,1.8,2.]:
        x=tri_evolve(mu,[.2,.25],ts);tracks[f'triangle_{mu}']=x
        if mu<np.pi/2:
            pred=0.;numeric=0.;phase=mu
        else:
            period=quad(lambda p:1/(mu-float(triangle(p))),-np.pi/2,3*np.pi/2,points=[np.pi/2],epsabs=1e-11)[0]
            numeric=2*np.pi/period;pred=np.pi/np.log((mu+np.pi/2)/(mu-np.pi/2));phase=None
        tri.append({'mu':mu,'stable_phase':phase,'predicted_drift':float(pred),'quadrature_drift':float(numeric)})
    As={
        'stable_node':np.diag([-2.,-1.]),'saddle':np.diag([1.,-1.]),
        'stable_spiral':np.array([[-.3,-1.],[1.,-.3]]),'center':np.array([[0.,2.],[-1.,0.]]),
        'degenerate_line':np.diag([0.,-1.]),'nonnormal_stable':np.array([[-1.,8.],[0.,-2.]])}
    times=np.linspace(0,5,501);linear=[]
    for name,A in As.items():
        maps=np.array([expm(A*t) for t in times]);gain=np.linalg.svd(maps,compute_uv=False)[:,0]**2
        sol=solve_ivp(lambda t,x:A@x,(0,5),[.3,.7],t_eval=times,rtol=1e-11,atol=1e-13,method='DOP853').y.T
        exact=np.einsum('tij,j->ti',maps,[.3,.7]);tracks[name+'_path']=sol;tracks[name+'_gain']=gain
        eig=np.linalg.eigvals(A)
        linear.append({'name':name,'matrix':A.tolist(),'trace':float(np.trace(A)),'determinant':float(np.linalg.det(A)),
                       'eigenvalues':[[float(x.real),float(x.imag)] for x in eig],
                       'classification':classify(A),'maximum_gain_0_to_5':float(gain.max()),'peak_time':float(times[np.argmax(gain)]),
                       'numerical_vs_exponential_max_error':float(abs(sol-exact).max())})
    center=tracks['center_path'];invariant=center[:,0]**2+2*center[:,1]**2
    np.savez_compressed(out/'dynamics.npz',linear_time=times,triangle_time=ts,**tracks)
    save_json(out/'dynamics.json',{'triangle':tri,'linear':linear,'center_max_invariant_drift':float(abs(invariant-invariant[0]).max()),
            'triangle_raw_lock_range':'|Delta| < A*pi/2; phi*=Delta/A on the stable rising branch',
            'triangle_normalized_lock_range':'For g=2f/pi, |Delta|<A; the wider raw range is an amplitude-normalization effect',
            'circuit_status':'Conditional overdamped identical-junction RSJ equations and algebraic tests only. No component-specific circuit experiment is claimed.'})
    print('Triangular-response and linear-system numerical controls complete',flush=True)


def supplement():
    from analyze import FIG,RESULTS,plt,COLORS,figsave,img,table,read
    r=read(RESULTS/'dynamics.json');z=np.load(RESULTS/'dynamics.npz')
    p=np.linspace(-np.pi,2*np.pi,600);fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    axs[0].plot(p,np.sin(p),label='sin φ');axs[0].plot(p,triangle(p),label='Raw triangle f');axs[0].plot(p,2/np.pi*triangle(p),'--',label='Equal-peak triangle 2f/π')
    axs[0].set(title='Response shape and normalization',xlabel='φ (rad)',ylabel='Coupling response');axs[0].legend(fontsize=8)
    mus=np.linspace(-.999,.999,300);axs[1].plot(mus,np.arcsin(mus),label='Sine');axs[1].plot(mus,np.pi/2*mus,label='Equal-peak triangle')
    axs[1].set(title='Locked phase at equal peak forcing',xlabel='Δ / A',ylabel='Stable phase (rad)');axs[1].legend(fontsize=8)
    for mu in [0.,.5,1.4,1.8,2.]:
        x=z[f'triangle_{mu}'];axs[2].semilogy(z['triangle_time'],abs(x[:,1]-x[:,0]),label=f'μ={mu}')
    axs[2].set(title='Raw triangular phase-kick response',xlabel='Normalized time',ylabel='Unwrapped phase separation');axs[2].legend(fontsize=8)
    for ax in axs:ax.grid(alpha=.2)
    figsave(fig,'triangular_response')
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    dd=np.linspace(-2,5,501);tt=np.linspace(-5,5,501);D,T=np.meshgrid(dd,tt)
    disc=T*T-4*D;types=np.where(D<0,0,np.where(disc>=0,np.where(T<0,1,2),np.where(T<0,3,4)))
    from matplotlib.colors import ListedColormap
    axs[0].imshow(types,origin='lower',extent=[-2,5,-5,5],aspect='auto',cmap=ListedColormap(['#dcc9ad','#c7e0d6','#f0d0cc','#b5d4e7','#e3d4eb']))
    a=np.linspace(-4.47,4.47,500);axs[0].plot(a*a/4,a,color='.4',lw=.8);axs[0].axhline(0,color='.4',lw=.8);axs[0].axvline(0,color='.4',lw=.8)
    for x,y,label in [(-1,0,'Saddles'),(1,-3.9,'Stable nodes'),(1,3.9,'Unstable nodes'),(3,-1.5,'Stable spirals'),(3,1.5,'Unstable spirals'),(3,0,'Centers')]:axs[0].text(x,y,label,ha='center',fontsize=8)
    axs[0].set(title='2 × 2 autonomous classification',xlabel='Determinant Δ',ylabel='Trace τ')
    for name,label in [('stable_node','Normal stable matrix'),('nonnormal_stable','Same eigenvalues; off-diagonal coupling')]:
        axs[1].plot(z['linear_time'],z[name+'_gain'],label=label)
    axs[1].axhline(1,color='.4',ls=':');axs[1].set(title='Stable eigenvalues can hide transient growth',xlabel='Time',ylabel='Optimal squared-norm gain');axs[1].legend(fontsize=8)
    for point in [[1,0],[0,1],[1,1],[-1,.5]]:
        A=np.array([[0.,2.],[-1.,0.]]);path=np.array([expm(A*t)@point for t in np.linspace(0,5,250)])
        axs[2].plot(path[:,0],path[:,1])
    axs[2].set(title='Undamped reciprocal coupling: center',xlabel='R',ylabel='J',aspect='equal');axs[2].grid(alpha=.2)
    figsave(fig,'linear_dynamics')
    # Reproduce the a-family qualitative portraits using an explicit illustrative model.
    fig,axs=plt.subplots(1,5,figsize=(13,2.6),layout='constrained');axis=np.linspace(-1.5,1.5,50);X,Y=np.meshgrid(axis,axis)
    for ax,a in zip(axs,[-2.,-1.,-.4,0.,.7]):
        ax.streamplot(axis,axis,a*X,-Y,density=.65,color=COLORS[0],linewidth=.7,arrowsize=.8)
        ax.set(title=f'a={a:g}',xlabel='x',ylabel='y',aspect='equal');ax.set_xticks([]);ax.set_yticks([])
        if a==0:ax.plot(axis,np.zeros_like(axis),color=COLORS[1],lw=2)
    fig.suptitle('Representative portrait family: ẋ = ax, ẏ = −y');figsave(fig,'phase_portraits')
    nm=next(x for x in r['linear'] if x['name']=='nonnormal_stable')
    parts=['<h2>7b. Added excerpts: response shape, stability and coupled systems</h2><h3>Triangular firefly response</h3><p>The second excerpt replaces the sine response by a periodic triangular wave: f(φ)=φ on [−π/2,π/2] and f(φ)=π−φ on [π/2,3π/2]. For the same phase convention, φ̇=Δ−Af(φ). The stable branch has φ*=Δ/A, provided |Δ|&lt;Aπ/2. Its small-perturbation recovery rate is A. At the endpoints the branches meet at a nonsmooth corner; the smooth sine model’s square-root critical slowing does not transfer unchanged.</p><p>The raw triangle has peak π/2, while sinφ has peak 1. Its larger locking range at equal A therefore partly reflects larger forcing. After equal-peak normalization g=2f/π, both locking intervals are |Δ|&lt;A; the equilibrium phase and relaxation rate still differ. Above the raw threshold, the analytic normalized mean drift is π / log[(μ+π/2)/(μ−π/2)] for μ&gt;π/2. Numerical quadrature and paired trajectories verify these controls.</p>',img('triangular_response','A response-shape null model: normalize forcing amplitude before claiming stronger entrainment.'),
        '<h3>Eigenvalues, nodes, saddles and the trace–determinant plane</h3><div class="eq">ẋ=Ax, &nbsp; λ²−τλ+Δ=0, &nbsp; τ=tr(A), Δ=det(A).<br>λ± = [τ ± √(τ²−4Δ)]/2.</div><p>For a real autonomous 2×2 system, Δ&lt;0 gives a saddle. With Δ&gt;0, negative trace gives decay and positive trace gives growth; real eigenvalues produce nodes and complex eigenvalues spirals. Trace zero with Δ&gt;0 gives a center for the linear system. Δ=0 or a repeated eigenvalue needs separate inspection: eigenvectors/Jordan structure determine a line of equilibria, a star node or a defective node. A nonlinear system with purely imaginary or zero eigenvalues cannot be classified from linearization alone.</p>',
        img('phase_portraits','These are explicit illustrative linear systems matching the supplied qualitative portrait family, not measured fluid fields.'),
        f'<p>The completed matrix control A=[[-1,8],[0,-2]] has eigenvalues −1 and −2 yet reaches optimal squared-norm gain <b>{nm["maximum_gain_0_to_5"]:.3f}</b> at t≈{nm["peak_time"]:.3f}. The diagonal matrix with the same eigenvalues only decays. The off-diagonal coupling makes the first matrix non-normal; asymptotic stability alone cannot rule out finite-time amplification. These numbers come from matrix exponentials verified against numerical ODE integration.</p>',
        img('linear_dynamics','Trace–determinant classification describes asymptotic eigenvalue behavior; finite-time growth requires the propagator or an energy estimate.'),
        '<p>The reciprocal example Ṙ=aJ, J̇=−bR with a,b&gt;0 has eigenvalues ±i√(ab) and conserved bR²+aJ². The simulated case a=2, b=1 produces elliptical center orbits. These equations are a mathematical oscillator example; the variable names in the excerpt do not make it a validated model of people.</p>',
        '<h3>Connection to fluid organization</h3><p>For neighboring passive trajectories, δẊ≈[∇u(X,t)]δX and d|δX|²/dt=2δXᵀSδX, where S is the symmetric velocity gradient. Finite-time deformation depends on the ordered history of these matrices; a single instantaneous eigenvalue need not predict it. This is one established reason that past geometry can help a partial-observation forecast. Material separation is still distinct from the global fluid perturbation-energy budget in Section 2.</p><p>The trace–determinant diagram is two-dimensional and autonomous. The pilot is three-dimensional and time-dependent. Incompressibility requires tr(∇u)=0, so arbitrary compressing/expanding phase-plane nodes cannot be assigned to local fluid-volume dynamics. Use the full 3D strain, deformation map and projected perturbation operator to test fluid claims.</p>',
        '<h3>Shared-load circuit: a conditional coupling model</h3><p>The circuit excerpt depicts two phase-labeled junction elements, individual shunts r, a shared load R and a bias current. If these are identical overdamped Josephson junctions with critical current I<sub>c</sub>, negligible capacitance/inductance and V<sub>i</sub>=(ℏ/2e)φ̇<sub>i</sub>, Kirchhoff’s law and I<sub>a</sub>=I<sub>c</sub>sinφ<sub>i</sub>+V<sub>i</sub>/r give the conditional dimensionless equations below. The phase–voltage relation is documented in <a href="https://www.nist.gov/document/voltage-metrology-superconductive-electronics-presentation">NIST’s Josephson-device presentation</a>.</p><div class="eq">(1+α)φ₁′ + αφ₂′ = i<sub>b</sub>−sinφ₁<br>αφ₁′ + (1+α)φ₂′ = i<sub>b</sub>−sinφ₂,<br>α=r/R, &nbsp; i<sub>b</sub>=I<sub>b</sub>/I<sub>c</sub>, &nbsp; τ=(2eI<sub>c</sub>r/ℏ)t.</div><p>This shows how an environmental load transmits coupling; R→∞ removes this shared-load term. A runnable RHS and symmetry/decoupling tests are included. No component-specific circuit simulation or synchronization conclusion is claimed because the attachment supplies no critical currents, capacitances or operating point. Coupled equations alone do not guarantee attracting phase synchronization.</p>',
        '<details><summary>Additional supplied reference images</summary>']
    for f in sorted((ROOT/'references').glob('codex-clipboard*.png')):parts.append(f'<img src="references/{f.name}" alt="User-supplied dynamical-systems reference excerpt">')
    parts.append('</details>')
    return '\n'.join(parts)


if __name__=='__main__':main()
