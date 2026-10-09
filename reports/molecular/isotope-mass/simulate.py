"""Confined explicit TIP3P isotope-mass experiment. Units: nm, ps, dalton, kJ/mol.
Run: python simulate.py --out data --workers 4
No bulk-water, nuclear-quantum or Navier-Stokes interpretation is implied.
"""
from pathlib import Path
import argparse, concurrent.futures, hashlib, json, os, platform, sys, time
import numpy as np
import openmm as mm
from openmm import app, unit

BASE = 20261009
H, D, O = 1.007947, 2.01410177812, 15.99943
R_GAS = .00831446261815324

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()),encoding='utf-8')

def system(n=9, heavy=-1, radius=None, mass_scale=None):
    top=app.Topology(); chain=top.addChain()
    for i in range(n):
        r=top.addResidue('HOH',chain,str(i)); a=[top.addAtom(name,el,r) for name,el in [('O',app.element.oxygen),('H1',app.element.hydrogen),('H2',app.element.hydrogen)]]
        top.addBond(a[0],a[1]);top.addBond(a[0],a[2])
    s=app.ForceField('tip3p.xml').createSystem(top,nonbondedMethod=app.NoCutoff,rigidWater=True,removeCMMotion=False)
    for i in range(n):
        s.setParticleMass(3*i,O)
        for j in [1,2]:s.setParticleMass(3*i+j,D if i==heavy else H)
    if mass_scale is not None and heavy>=0:
        for j,m in enumerate([O,H,H]):s.setParticleMass(3*heavy+j,m*mass_scale)
    radius = radius if radius is not None else .4*(n/9)**(1/3)
    w=mm.CustomExternalForce('0.5*k*max(0,sqrt(x*x+y*y+z*z)-R)^2')
    w.addGlobalParameter('k',10000.);w.addGlobalParameter('R',radius)
    for i in range(n):w.addParticle(3*i,[])
    s.addForce(w)
    return s,top

def masses(s):return np.array([s.getParticleMass(i).value_in_unit(unit.dalton) for i in range(s.getNumParticles())])

def random_positions(n, config, radius):
    rng=np.random.default_rng(BASE+1000*n+config); pts=[]
    for attempt in range(1000000):
        q=rng.uniform(-radius,radius,3)
        if np.linalg.norm(q)<radius and all(np.linalg.norm(q-p)>.235 for p in pts):pts.append(q)
        if len(pts)==n:break
    else:raise RuntimeError('packing failed')
    # Site index denotes radial rank of the randomized INITIAL packing, not a fixed spatial seat.
    pts=np.array(sorted(pts,key=np.linalg.norm)); out=[]
    water=np.array([[0,0,0],[.09572,0,0],[.09572*np.cos(np.deg2rad(104.52)),.09572*np.sin(np.deg2rad(104.52)),0]])
    for p in pts:
        q,r=np.linalg.qr(rng.normal(size=(3,3)));q=q@np.diag(np.sign(np.diag(r)))
        if np.linalg.det(q)<0:q[:,0]*=-1
        out.extend(water@q+p)
    return np.array(out)

def context(s,dt=.001,gamma=1,seed=1,nve=False):
    integ=mm.VerletIntegrator(dt) if nve else mm.LangevinMiddleIntegrator(300,gamma,dt)
    integ.setConstraintTolerance(1e-8)
    if not nve:integ.setRandomNumberSeed(seed)
    return mm.Context(s,integ,mm.Platform.getPlatformByName('Reference')),integ

def get_state(ctx):
    st=ctx.getState(getPositions=True,getVelocities=True,getEnergy=True)
    return (st.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
            st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond),
            [st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)])

def response(s,q,v,dt=.001,impulse=.05,sources=None):
    """Central +/- equal-total-momentum molecular kicks, 3 Cartesian axes, NVE.
    Save every branch position so the lagged intervention assay can be reanalyzed.
    No thermostat or momentum removal in the response window.
    """
    n=len(q)//3; sources=list(range(n)) if sources is None else list(sources)
    lags=np.array([0,.02,.1,.5]); steps=np.rint(lags/dt).astype(int)
    ctx,integ=context(s,dt,nve=True); m=masses(s).reshape(n,3).sum(axis=1)
    raw=np.zeros((len(sources),3,2,len(lags),n,3))
    base=np.zeros((len(lags),n,3)); baseline_energy=[]
    ctx.setPositions(q);ctx.setVelocities(v);ctx.applyVelocityConstraints(1e-8)
    for l in range(len(lags)):
        if l:integ.step(int(steps[l]-steps[l-1]))
        a,_,e=get_state(ctx);base[l]=a[::3];baseline_energy.append(sum(e))
    for si,src in enumerate(sources):
        for axis in range(3):
            for sign_id,sign in enumerate([-1,1]):
                kick=v.copy();kick[3*src:3*src+3,axis]+=sign*impulse/m[src]
                ctx.setPositions(q);ctx.setVelocities(kick);ctx.applyVelocityConstraints(1e-8)
                for l in range(len(lags)):
                    if l:integ.step(int(steps[l]-steps[l-1]))
                    raw[si,axis,sign_id,l]=ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)[::3]
    del ctx,integ
    return raw,base,np.array(baseline_energy),lags

def tasks():
    specs=[]
    # Each (configuration, velocity seed) is a paired block; configurations are bootstrap clusters.
    suites=[('main',9,16,2,50,100,.001,1,.4,list(range(9)),True),
            ('long',9,8,2,100,300,.001,1,.4,[0,8],False),
            ('size27',27,8,2,100,100,.001,1,.4*3**(1/3),[0,26],False),
            ('wall',9,8,2,50,100,.001,1,.5,[0,8],False),
            ('halfdt',9,4,2,50,100,.0005,1,.4,[0,8],True),
            ('friction',9,4,2,50,100,.001,.1,.4,[0,8],False)]
    for suite,n,nconfig,nvel,eq,prod,dt,gamma,radius,sites,resp in suites:
        for c in range(nconfig):
            for v in range(nvel):
                for heavy in [-1]+sites:
                    specs.append(dict(suite=suite,n=n,config=c,velocity=v,heavy=heavy,equil_ps=eq,production_ps=prod,dt_ps=dt,gamma_ps=gamma,radius_nm=radius,response=resp,label=False))
                if suite=='main':
                    specs.append(dict(specs[-len(sites)-1],label=True))
    return specs

def run_one(args):
    task,out=args; out=Path(out)
    stem=f"{task['suite']}_c{task['config']:02d}_v{task['velocity']}_"+('label' if task['label'] else 'A' if task['heavy']<0 else f"D{task['heavy']:02d}")
    if (out/(stem+'.json')).exists() and (out/(stem+'.npz')).exists():
        previous=json.loads((out/(stem+'.json')).read_text())
        if any(previous.get(k)!=v for k,v in task.items()):
            raise ValueError(stem+' has different saved settings; choose an empty output directory')
        return stem+' cached'
    tic=time.time();n=task['n'];s,top=system(n,task['heavy'],task['radius_nm']);m=masses(s)
    seed=BASE+task['config']*100+task['velocity']*17
    initial=random_positions(n,task['config'],task['radius_nm'])
    ctx,integ=context(s,task['dt_ps'],task['gamma_ps'],seed)
    ctx.setPositions(initial);mm.LocalEnergyMinimizer.minimize(ctx,1.,2000)
    minimized=get_state(ctx)[0]
    # Matched Gaussian variates across arms; mass-dependent Maxwell velocities.
    rng=np.random.default_rng(seed+7)
    vel=rng.normal(size=(3*n,3))*np.sqrt(R_GAS*300/m)[:,None]
    ctx.setVelocities(vel);ctx.applyVelocityConstraints(1e-8)
    eq_trace=[]
    for _ in range(round(task['equil_ps'])):
        integ.step(round(1/task['dt_ps']));_,_,e=get_state(ctx);eq_trace.append(e)
    sample_ps=.05; nf=round(task['production_ps']/sample_ps)+1
    coords=np.empty((nf,3*n,3),dtype=np.float32);energies=np.empty((nf,2));snap_q=[];snap_v=[]
    snap_frames={round((nf-1)*f) for f in [.25,.5,.75]}
    for f in range(nf):
        if f:integ.step(round(sample_ps/task['dt_ps']))
        q,v,e=get_state(ctx);coords[f]=q;energies[f]=e
        if f in snap_frames:snap_q.append(q);snap_v.append(v)
    perm=np.random.default_rng(seed+99).permutation(n) if task['label'] else np.arange(n)
    # Labels are metadata only. Reordering stored molecule rows changes no forces/masses/noise.
    coords=coords.reshape(nf,n,3,3)[:,perm]
    raw=dict(positions_nm=coords,energies_kjmol=energies,equil_energies_kjmol=eq_trace,
             initial_nm=initial.reshape(n,3,3),minimized_nm=minimized.reshape(n,3,3),
             snapshot_q_nm=snap_q,snapshot_v_nmps=snap_v,permutation=perm,masses_da=m.reshape(n,3)[perm])
    if task['response'] and not task['label']:
        rr=[];bb=[];ee=[]
        for q,v in zip(snap_q,snap_v):
            a,b,e,lags=response(s,q,v,task['dt_ps']);rr.append(a);bb.append(b);ee.append(e)
        raw.update(response_positions_nm=np.array(rr),response_baseline_nm=np.array(bb),response_energy_kjmol=np.array(ee),lags_ps=lags,impulse_da_nmps=.05)
    if not np.isfinite(coords).all() or not np.isfinite(energies).all():raise ValueError(stem+' nonfinite')
    np.savez_compressed(out/(stem+'.npz'),**raw)
    task=dict(task,seed=seed,sample_ps=sample_ps,frames=nf,seconds=time.time()-tic,platform='Reference',temperature_K=300,openmm=mm.__version__)
    dump(out/(stem+'.json'),task)
    del ctx,integ
    return stem+f" {task['seconds']:.1f}s"

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('data'));p.add_argument('--workers',type=int,default=4);p.add_argument('--suite',default='all');p.add_argument('--limit',type=int,default=0);args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True); jobs=[t for t in tasks() if args.suite=='all' or t['suite']==args.suite]
    if args.limit:jobs=jobs[:args.limit]
    dump(args.out/'run_environment.json',dict(python=sys.version,os=platform.platform(),openmm=mm.__version__,numpy=np.__version__,argv=sys.argv,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    for n in [9,27]:
        for heavy in [-1,0]:
            s,top=system(n,heavy);(args.out/f'system_n{n}_h{heavy}.xml').write_text(mm.XmlSerializer.serialize(s),encoding='utf-8')
    dump(args.out/'planned_tasks.json',jobs)
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        for i,msg in enumerate(pool.map(run_one,[(t,str(args.out)) for t in jobs])):print(f'{i+1}/{len(jobs)} {msg}',flush=True)

if __name__=='__main__':main()
