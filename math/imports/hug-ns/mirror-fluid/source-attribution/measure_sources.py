"""Controlled source additions to the supplied AB field; unchanged fluid solver.

At t=.5, use inclusion-exclusion on the signed fluid observable S.
S0 = S(AB); delta_C = S(AB+aC)-S0; delta_D = S(AB+bD)-S0;
interaction = S(AB+aC+bD)-S(AB+aC)-S(AB+bD)+S0.
This defines an endpoint interaction contrast, not a unique decomposition of
the evolving nonlinear velocity, and not fitted terms in the proposed q ODE.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source'))
import clay_hug as m
import fluid_sources as original
m.NU=.01
N=33
END=.5
CASES={'half_c':(.5,0.),'half_d':(0.,.5),'seven_c':(.7,0.),'three_d':(0.,.3),
       'equal':(.5,.5),'unequal':(.7,.3)}
JOBS=[{'id':k+'-base','case':k,'half':False} for k in ('half_c','half_d','seven_c','three_d')]+[
     {'id':k+'-half','case':k,'half':True} for k in CASES]

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8');tmp.replace(p)
def hashes():
    return {str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
       (HERE/'measure_sources.py',HERE/'source/clay_hug.py',HERE/'source/fluid_sources.py')}
def build():
    x,y,z,kx,ky,kz,k2,keep,dx=m.grids(N);ops=(kx,ky,kz,k2,keep,dx)
    psi=np.exp(-12*((x-.9)**2+y*y))*(x>0)
    ab=[m.deriv(psi,ky,keep),-m.deriv(psi,kx,keep),np.zeros_like(psi)]
    ab=[(a+b)/2 for a,b in zip(ab,original.mirror_of(ab))]
    ab=m.project(ab,kx,ky,kz,k2,keep)
    ab=[(a+b)/2 for a,b in zip(ab,original.mirror_of(ab))]
    blob=np.exp(-8*((x-.55)**2+(y-.4)**2+z*z))
    c=m.project([np.zeros_like(blob),.25*blob,np.zeros_like(blob)],kx,ky,kz,k2,keep)
    d=original.mirror_of(c)
    phi=[(a-b)/2 for a,b in zip(c,d)]
    phi=[a/np.sqrt(original.dot(phi,phi)) for a in phi]
    return ab,c,d,phi,ops
def field(ab,c,d,weights,ops):
    a,b=weights
    return m.project([x+a*y+b*z for x,y,z in zip(ab,c,d)],*ops[:-1])
def signed(u,phi):return original.dot([a-b for a,b in zip(u,original.mirror_of(u))],phi)
def observe(u,t,phi,ops):
    kx,ky,kz,k2,keep,dx=ops
    uh=np.array([m.spec(v,keep) for v in u])
    oh=np.array([1j*(ky*uh[2]-kz*uh[1]),1j*(kz*uh[0]-kx*uh[2]),1j*(kx*uh[1]-ky*uh[0])])
    modes=np.fft.fftfreq(N)*N
    high=((abs(modes[:,None,None])>=.8*(N//3))|(abs(modes[None,:,None])>=.8*(N//3))|(abs(modes[None,None,:])>=.8*(N//3)))&keep
    return {'t':t,'S':signed(u,phi),'E':original.em(u),'W':original.wmax(u,ops),
      'energy':m.energy(u),'divergence':m.div_max(u,kx,ky,kz,keep),
      'high_band_energy_fraction':float(np.sum(abs(uh[:,high])**2)/np.sum(abs(uh)**2)),
      'high_band_enstrophy_fraction':float(np.sum(abs(oh[:,high])**2)/np.sum(abs(oh)**2))}
def setup():
    ab,c,d,phi,ops=build();ref=read(HERE/'reference-fluid-sources.json');checks=[]
    for label,weights in [('ab',(0,0)),('ab_cd',(.5,.5)),('ab_c',(1,0)),('ab_uneven',(.7,.3)),('ab_c_mirror',(0,1))]:
        u=field(ab,c,d,weights,ops);row=observe(u,0.,phi,ops)
        errors={k:abs(row[k]-ref[label][0][alias]) for k,alias in [('S','D'),('E','E'),('W','W')]}
        assert max(errors.values())<1e-10,(label,errors)
        checks.append({'case':label,'initial_errors':errors})
    save(HERE/'protocol.json',{'purpose':'Measure individual source additions and their interaction at the existing endpoint',
      'n':N,'end':END,'nu':m.NU,'domain_side':m.L,'force':0,'method':'Unchanged supplied Heun solver and projection',
      'starting_field':'Exact supplied fluid_sources.py AB, C, and D=mirror(C)',
      'template':'Same normalized antisymmetric part of C as fluid_sources.py',
      'jobs':JOBS,'weights':CASES,'reference_sha256':hashlib.sha256((HERE/'reference-fluid-sources.json').read_bytes()).hexdigest(),
      'reuse':'Five completed source trajectories are retained. Only missing single-source amplitudes and timestep controls are evolved.',
      'observables':{'S':'<u-M[u],phi>; called D in the supplied JSON','E':'norm(u-M[u])/norm(u)',
        'interaction':'S(AB+aC+bD)-S(AB+aC)-S(AB+bD)+S(AB), evaluated at t=0.5'},
      'limits':['Initial weights are controlled amplitudes, not fitted q couplings.',
        'A nonzero interaction contrast detects nonadditivity, not automatically cross-pair competition.',
        'E is a nonlinear unsigned norm and is not decomposed by addition.',
        'AB is held mirror-symmetric; these controls do not identify an independent AB-bias response.',
        'D is the exact mirror of C, not an independently shaped source. No grid refinement is included in this bounded attribution test.'],
      'source_hashes':hashes(),'initial_checks':checks})
def run(job):
    folder=HERE/'runs'/job['id'];p=folder/'result.json'
    if p.exists() and read(p)['status']=='complete':return read(p)
    ab,c,d,phi,ops=build();u=field(ab,c,d,CASES[job['case']],ops)
    speed=float(np.sqrt(sum(v*v for v in u)).max())
    base_steps=int(np.ceil(END/(.05*ops[-1]/max(speed,1e-6))))
    factor=2 if job['half'] else 1;steps=base_steps*factor;dt=END/steps
    every=max(1,base_steps//5)*factor;series=[];start=0;max_cfl=0.
    if p.exists():
        old=read(p);assert old['source_hashes']==hashes() and old['dt']==dt
        series=old['series'];start=old['last_step']+1;max_cfl=old['max_speed_cfl']
        u=list(np.load(folder/old['field']))
    for step in range(start,steps+1):
        if step:u=m.advance(u,dt,ops)
        assert np.isfinite(u).all(),(job,'nonfinite',step)
        max_cfl=max(max_cfl,float(np.sqrt(sum(v*v for v in u)).max())*dt/ops[-1])
        assert max_cfl<.5,(job,'CFL',max_cfl)
        if step%every==0 or step==steps:
            row=observe(u,step*dt,phi,ops)
            assert row['divergence']<1e-10,(job,'divergence',row)
            series.append(row);folder.mkdir(parents=True,exist_ok=True)
            name=f'field-step{step:05d}.npy';np.save(folder/name,np.array(u))
            result={'status':'complete' if step==steps else 'running','job':job,'n':N,'weights':CASES[job['case']],
              'dt':dt,'steps':steps,'last_step':step,'field':name,'max_speed_cfl':max_cfl,'source_hashes':hashes(),'series':series}
            save(p,result)
            print(f"{job['id']} t={row['t']:.5f} S={row['S']:.9g} E={row['E']:.9g}",flush=True)
    return result
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    setup()
    state={'status':'running','pid':os.getpid(),'planned':len(JOBS),'completed':[],'failures':{},'started_epoch':time.time()}
    save(HERE/'task-state.json',state)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run,j):j['id'] for j in JOBS}
        for f in as_completed(futures):
            name=futures[f]
            try:
                f.result();state['completed'].append(name)
            except Exception as e:state['failures'][name]=repr(e)
            save(HERE/'task-state.json',state)
    state['status']='complete' if not state['failures'] else 'failed';save(HERE/'task-state.json',state)
    print(json.dumps(state),flush=True)
    if state['failures']:raise SystemExit(1)
if __name__=='__main__':main()
