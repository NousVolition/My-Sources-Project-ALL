"""Fluid data for a measured-asymmetry closure test. Does not force the fluid.

Fixed six-template projection defines candidate source features. Fitting uses
only training trajectories. q_fluid is measured; q_record is an imposed linear
diagnostic and cannot independently establish nonlinear feedback.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse, hashlib, json, os, sys, time
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source'))
import clay_hug as m
import fluid_sources as original
m.NU=.01
END=.5
TIMES=np.linspace(0,END,51)
TRAIN=[('zero',0,0,0),('ab_plus',.01,0,0),('ab_minus',-.01,0,0),
       ('c_only',0,1,0),('d_only',0,0,1),('equal',0,.5,.5),
       ('mix_a',.005,.7,.3),('mix_b',-.005,.3,.7),
       ('mix_c',.01,1,.5),('mix_d',-.01,.5,1)]
TEST=[('unseen_a',.003,.8,.2,0),('unseen_b',-.007,.2,.8,0),
      ('unseen_shape_a',.006,.6,.4,1),('unseen_shape_b',-.004,.4,.6,2)]
CASES={name:dict(bias=a,c=c,d=d,shape=0,split='train') for name,a,c,d in TRAIN}
CASES.update({name:dict(bias=a,c=c,d=d,shape=s,split='test') for name,a,c,d,s in TEST})
JOBS=[dict(id=name+'-n33',case=name,n=33,half=False,mirror=False) for name in CASES]
JOBS += [dict(id=name+suffix,case=name,n=n,half=half,mirror=mirror)
         for suffix,n,half,mirror in [('-mirror',33,False,True),('-half',33,True,False),('-n49',49,False,False)]
         for name,*_ in TEST]

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for attempt in range(20):
        try:tmp.replace(p);return
        except PermissionError:
            if attempt==19:raise
            time.sleep(.05)
def hashes():return {str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
    for p in [HERE/'study.py',HERE/'source/clay_hug.py',HERE/'source/fluid_sources.py']}
def mirror(u):return np.asarray(original.mirror_of(u))
def dot(a,b):return float(np.sum(a*b)/np.prod(a.shape[1:]))
def norm(a):return np.sqrt(dot(a,a))
def unit(a):return a/norm(a)
def parity(u,p):return (u+p*mirror(u))/2

def build(n):
    x,y,z,kx,ky,kz,k2,keep,dx=m.grids(n);ops=(kx,ky,kz,k2,keep,dx)
    def project(a):return np.asarray(m.project(a,*ops[:-1]))
    psi=np.exp(-12*((x-.9)**2+y*y))*(x>0)
    ab=parity(project(parity(np.asarray([m.deriv(psi,ky,keep),-m.deriv(psi,kx,keep),np.zeros_like(psi)]),1)),1)
    blob=np.exp(-8*((x-.55)**2+(y-.4)**2+z*z))
    c=project([np.zeros_like(blob),.25*blob,np.zeros_like(blob)])
    blob=np.exp(-8*((x-.9)**2+(y-.35)**2+z*z))*(x>0)
    bias=unit(parity(project([np.zeros_like(blob),blob,np.zeros_like(blob)]),-1))
    ds=[]
    for xc,yc,zc,width in [(.7,.15,.2,6.),(.8,.1,.25,5.),(.62,.25,-.12,7.)]:
        blob=np.exp(-width*((x-xc)**2+(y-yc)**2+(z-zc)**2))
        ds.append(mirror(project([np.zeros_like(blob),.25*blob,.1*blob])))
    scale=norm(ab);phi=unit(parity(c,-1))
    basis=np.asarray([unit(ab),bias,unit(parity(c,1)),phi,unit(parity(ds[0],1)),unit(parity(ds[0],-1))])
    flat=basis.reshape(6,-1);gram=flat@flat.T/n**3
    inverse=np.linalg.inv(gram)
    # B_ij = projection of -P(e_i . grad e_j) onto the fixed signed template.
    # phi is solenoidal and retained by the filter, so P drops from this inner product.
    tensor=np.zeros((6,6))
    for j in range(6):
        grad=np.asarray([[m.deriv(basis[j,k],kk,keep) for kk in (kx,ky,kz)] for k in range(3)])
        for i in range(6):
            adv=np.einsum('bxyz,abxyz->axyz',basis[i],grad,optimize=True)
            tensor[i,j]=-dot(adv,phi)/scale
    blocks=np.zeros((4,6,6))
    for i in range(6):
        for j in range(6):
            a,b=i//2,j//2
            block=3 if {a,b}=={1,2} else max(a,b)
            blocks[block,i,j]=tensor[i,j]
    signs=np.array([1,-1,1,-1,1,-1])
    even_mask=signs[:,None]*signs[None,:]==1
    assert np.max(abs(tensor[even_mask]))<1e-10
    # Direct independent check against the supplied full RHS at a fixed coefficient vector.
    coef=np.array([1.,.02,.03,-.01,.04,.015]);v=np.einsum('i,icxyz->cxyz',coef,basis)
    lap=np.asarray([np.fft.ifftn(-k2*m.spec(a,keep)).real for a in v])
    direct=dot(np.asarray(m.rhs(v,*ops[:-1]))-m.NU*lap,phi)/scale
    tensor_error=abs(direct-float(coef@tensor@coef))
    assert tensor_error<1e-10,(n,tensor_error)
    return dict(ab=ab,c=c,ds=ds,bias=bias,phi=phi,scale=scale,basis=basis,flat=flat,
                gram=gram,inverse=inverse,blocks=blocks,ops=ops,tensor_error=tensor_error)

def features(u,b):
    coeff=b['inverse']@(b['flat']@u.reshape(-1)/u.shape[1]**3)
    ds=np.einsum('i,kij,j->k',coeff,b['blocks'],coeff)
    residual=u-np.einsum('i,icxyz->cxyz',coeff,b['basis'])
    return ds,coeff,norm(residual)/norm(u)

def observe(u,t,qrec,b):
    ops=b['ops'];kx,ky,kz,k2,keep,dx=ops;n=u.shape[1]
    ds,coeff,residual=features(u,b)
    q=dot(u,b['phi'])/b['scale']
    rhs=np.asarray(m.rhs(u,*ops[:-1]))
    qdot=dot(rhs,b['phi'])/b['scale']
    uh=np.asarray([m.spec(v,keep) for v in u]);oh=np.asarray([1j*(ky*uh[2]-kz*uh[1]),1j*(kz*uh[0]-kx*uh[2]),1j*(kx*uh[1]-ky*uh[0])])
    modes=np.fft.fftfreq(n)*n
    high=((abs(modes[:,None,None])>=.8*(n//3))|(abs(modes[None,:,None])>=.8*(n//3))|(abs(modes[None,None,:])>=.8*(n//3)))&keep
    vorticity=np.asarray([np.fft.ifftn(v).real for v in oh])
    return dict(t=float(t),q_fluid=q,qdot_fluid=qdot,q_record=qrec,
        S=2*b['scale']*q,E=original.em(u),source_features=ds.tolist(),template_coefficients=coeff.tolist(),
        template_residual_fraction=residual,energy=m.energy(u),W=float(np.sqrt(np.sum(vorticity*vorticity,axis=0)).max()),
        divergence=m.div_max(u,kx,ky,kz,keep),
        high_band_energy_fraction=float(np.sum(abs(uh[:,high])**2)/np.sum(abs(uh)**2)),
        high_band_enstrophy_fraction=float(np.sum(abs(oh[:,high])**2)/np.sum(abs(oh)**2)))

def protocol():
    return dict(status='specified before fitting or viewing held-out trajectories',end=END,output_dt=.01,
        n_base=33,n_control=49,nu=.01,domain=6,force=0,method='Unchanged supplied Heun and projection',
        cases=CASES,jobs=JOBS,source_hashes=hashes(),
        initial='Supplied symmetric AB and C; independent signed AB template; D shifted, broadened and given a z component. D shape variants in held-out starts.',
        mirror='Reflect the entire initial vector field, with component signs; do not only flip a parameter.',
        coordinates={'q_fluid':'<u,phi>/norm(AB0) = S/(2 norm(AB0)); independently measured signed coordinate; not historical q_record',
          'q_record':'q_next=q+dt*(-0.2*q+0.8*S_next), q(0)=0. Imposed diagnostic control; not evidence for fluid feedback.',
          'E':'norm(u-M[u])/norm(u); unsigned, not inferred from q'},
        source_definition={'basis':['AB symmetric','independent AB antisymmetric bias','C symmetric','C antisymmetric (phi)','nominal independent D symmetric','nominal independent D antisymmetric'],
          'coefficients':'Least-squares spatial projection onto these fixed six templates; same basis for all cases at a grid.',
          'D_AB':'Signed projected nonlinear interactions within the two AB templates.',
          'D_C':'C self interactions and both ordered AB-C interactions.',
          'D_D':'D self interactions and both ordered AB-D interactions.',
          'D_cross':'Both ordered C-D interactions. This is an operational interaction feature, not an established competition-between-pairings mechanism.',
          'not_exact_closure':'Remaining spatial modes and full viscous projection are not supplied to the regressor; target full-fluid derivative is independently computed.'},
        fitting={'models':['decay_only','AB_only','AB_C','AB_C_D','full_linear','full_cubic'],
          'rule':'dq/dt=-gamma*q+sum(w_i*D_i)-g*q^3; gamma>=0,g>=0; w_i unrestricted',
          'ridge_candidates':[0.,1e-6,.001,.1], 'selection':'leave-one-training-run-out derivative error; training column scaling only',
          'unseen':'Four entire starts withheld; mirrors, timestep and grid controls excluded from fitting.',
          'conditional_prediction':'Integrate from initial q, supplying measured source features from held-out run, never future q. This needs the measured input trajectory.',
          'forecast':'At t=.1 retain observed q and source features, freeze inputs, predict to .5 without future observations; explicitly a persistence assumption.',
          'cubic_support':'Report full vs linear unseen errors, stability of coefficients across training omissions and refinements. Do not equate a fitted positive g with spontaneous branch selection.'},
        limits=['Candidate six-template definitions implement a measurable version of the proposed source terms; not unique causal source attribution.',
          '33 and49 are a grid sensitivity check, not a full spatial convergence ladder.',
          'With gamma>=0,g>=0 and constant source input, the scalar drift is monotone; it cannot have two isolated stable branches.',
          'No coefficient fit or selected prediction may alter the Navier-Stokes evolution.'])

def run(job):
    folder=HERE/'runs'/job['id'];p=folder/'result.json';current_hash=hashes()
    if p.exists() and read(p)['status']=='complete':
        old=read(p);assert old['source_hashes']==current_hash;return job['id']
    b=build(job['n']);case=CASES[job['case']];ops=b['ops']
    u=b['ab']+case['bias']*b['scale']*b['bias']+case['c']*b['c']+case['d']*b['ds'][case['shape']]
    u=np.asarray(m.project(u,*ops[:-1]))
    if job['mirror']:u=mirror(u)
    speed=float(np.sqrt(np.sum(u*u,axis=0)).max())
    per=int(np.ceil(.01/(.04*ops[-1]/speed)))*(2 if job['half'] else 1)
    dt=.01/per;series=[];qrec=0.;first=0;max_cfl=0.;integral=0.;previousW=None
    if p.exists():
        old=read(p);assert old['source_hashes']==current_hash and old['dt']==dt
        series=old['series'];u=np.load(folder/old['field']);qrec=series[-1]['q_record']
        first=old['last_output']+1;max_cfl=old['max_speed_cfl'];integral=old['sampled_I'];previousW=series[-1]['W']
    for index in range(first,len(TIMES)):
        if index:
            for _ in range(per):
                u=np.asarray(m.advance(u,dt,ops))
                assert np.isfinite(u).all(),(job['id'],'nonfinite')
                S=2*dot(u,b['phi']);qrec+=dt*(-.2*qrec+.8*S)
                max_cfl=max(max_cfl,float(np.sqrt(np.sum(u*u,axis=0)).max())*dt/ops[-1])
                assert max_cfl<.5,(job['id'],'CFL',max_cfl)
        row=observe(u,TIMES[index],qrec,b);assert row['divergence']*ops[-1]<1e-10
        if previousW is not None:integral+=.005*(row['W']+previousW)
        previousW=row['W'];row['sampled_I']=integral;series.append(row)
        folder.mkdir(parents=True,exist_ok=True);name=f'field-{index:03d}.npy'
        np.save(folder/name,u)
        result=dict(status='complete' if index==50 else 'running',job=job,case=case,dt=dt,steps_per_output=per,
          last_output=index,field=name,max_speed_cfl=max_cfl,sampled_I=integral,source_hashes=current_hash,
          scale=b['scale'],gram_condition=float(np.linalg.cond(b['gram'])),tensor_check_error=b['tensor_error'],series=series)
        save(p,result)
        if index%10==0:print(f"{job['id']} t={row['t']:.2f} q={row['q_fluid']:.8g} E={row['E']:.8g}",flush=True)
    return job['id']

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);parser.add_argument('--audit',action='store_true');args=parser.parse_args()
    if args.audit:
        b=build(33);u=b['ab']+.003*b['scale']*b['bias']+.8*b['c']+.2*b['ds'][1]
        a=observe(u,0,0,b);r=observe(mirror(u),0,0,b)
        assert abs(a['q_fluid']+r['q_fluid'])<1e-12
        assert abs(a['E']-r['E'])<1e-12
        assert max(abs(np.asarray(a['source_features'])+r['source_features']))<1e-10
        out=dict(tensor_error=b['tensor_error'],gram_condition=float(np.linalg.cond(b['gram'])),
            reflected_q_error=abs(a['q_fluid']+r['q_fluid']),reflected_E_error=abs(a['E']-r['E']),jobs=len(JOBS))
        save(HERE/'implementation-check.json',out);print(json.dumps(out));return
    lock=HERE/'controller.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    try:
        expected=protocol()
        if (HERE/'protocol.json').exists():assert read(HERE/'protocol.json')==expected
        else:save(HERE/'protocol.json',expected)
        state=dict(status='running',pid=os.getpid(),planned=len(JOBS),completed=[],failures={},started_epoch=time.time())
        save(HERE/'task-state.json',state)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures={pool.submit(run,j):j['id'] for j in JOBS}
            for f in as_completed(futures):
                name=futures[f]
                try:f.result();state['completed'].append(name)
                except Exception as e:state['failures'][name]=repr(e)
                save(HERE/'task-state.json',state)
        state['status']='complete' if not state['failures'] else 'failed';save(HERE/'task-state.json',state)
    finally:lock.unlink(missing_ok=True)
    if state['failures']:raise SystemExit(1)
if __name__=='__main__':main()
