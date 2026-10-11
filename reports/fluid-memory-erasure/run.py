"""Frozen-design runner. Raw recordings are resumed only after source/data checks."""
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np
from solver import Solver, random_field, curl_gaussian, interventions

ROOT = Path(__file__).resolve().parent
P = json.loads((ROOT/'protocol.json').read_text())
DIAG_NAMES = ['energy','energy_budget_relative','divergence_rms','tail_enstrophy','tangent_volume_error','gradient_trace','max_cfl']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def marker_cloud(seed):
    rng = np.random.default_rng(seed+170000)
    centers = np.array(P['probe_center'])+rng.uniform(-.7,.7,(P['centers'],3))
    offsets = np.concatenate([np.eye(3),-np.eye(3)])*P['neighbor_radius']
    clouds = []
    for point in centers:
        q,r = np.linalg.qr(rng.normal(size=(3,3)))
        if np.linalg.det(q)<0:
            q[:,0]*=-1
        clouds.append(point+offsets@q.T)
    return np.concatenate([centers,np.concatenate(clouds)])


def filename(seed,n,dt,nu,amplitude,end,backend):
    return f's{seed}-n{n}-dt{dt:g}-nu{nu:g}-p{amplitude:g}-T{end:g}-{backend}'


def run_one(seed,n=None,dt=None,nu=None,amplitude=None,end=None,backend='gpu',out=None):
    n = P['n'] if n is None else n
    dt = P['dt'] if dt is None else dt
    nu = P['nu'] if nu is None else nu
    amplitude = P['probe_velocity_rms'] if amplitude is None else amplitude
    end = max(P['horizons']) if end is None else end
    out = ROOT/'data' if out is None else Path(out)
    out.mkdir(parents=True,exist_ok=True)
    stem = filename(seed,n,dt,nu,amplitude,end,backend)
    meta_path,raw_path = out/(stem+'.json'),out/(stem+'.npz')
    sources = {f:digest(ROOT/f) for f in ['solver.py','run.py','protocol.json']}
    config = dict(seed=seed,n=n,dt=dt,nu=nu,amplitude=amplitude,end=end,backend=backend)
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if meta['sources']!=sources or meta['config']!=config or digest(raw_path)!=meta['raw_sha256']:
            raise RuntimeError('Cache fingerprint mismatch; use a fresh data folder')
        print(stem+' verified cache',flush=True)
        return meta
    tic = time.time()
    s = Solver(n,nu,backend)
    xp = s.xp
    initial = random_field(s,seed,P['initial_velocity_rms'])
    writer = curl_gaussian(s,P['probe_center'],.45,P['write_band_max'],high_only=True)
    probe = curl_gaussian(s,P['probe_center'],P['probe_width'],P['write_band_max'])*amplitude
    x0 = marker_cloud(seed)
    h = xp.stack([initial,initial])
    x = xp.asarray(np.stack([x0,x0]))
    tangent = xp.asarray(np.tile(np.eye(3),(2,P['centers'],1,1)))
    force = xp.stack([writer,-writer])*P['write_acceleration_rms']
    prep_energy0 = .5*s.inner(h,h)
    budget = xp.zeros(2)
    cflmax = xp.zeros(2)
    prept,prepx,prepu,prepg,preph,prepdiag = [],[],[],[],[],[]
    t = 0.
    for target in np.arange(0,P['prepare_end']+1e-9,P['save_dt']):
        while t < target-1e-12:
            step = min(dt,target-t)
            if t<P['write_stop']-1e-12:
                step = min(step,P['write_stop']-t)
                ff = force
            else:
                ff = xp.zeros_like(h)
            h,x,tangent,db,cfl = s.step(h,x,tangent,step,ff)
            budget+=db
            cflmax=xp.maximum(cflmax,cfl)
            t+=step
        u,g=s.sample(h,x,P['centers'])
        prept.append(float(target));prepx.append(s.cpu(x));prepu.append(s.cpu(u));prepg.append(s.cpu(g));preph.append(s.pack(h))
        prepdiag.append(s.diagnostics(h,x,tangent,budget,prep_energy0,cflmax))
    prepared = h.copy()
    variants,checks = interventions(s,prepared,seed,P['sham_high_band_translation'],P['phase_scramble_replicates'])
    variants = {'native':prepared,**variants}
    states,positions,labels = [],[],[]
    for arm,state in variants.items():
        for history in range(2):
            for kicked in range(2):
                states.append(state[history]+kicked*probe)
                positions.append(x[history] if arm=='native' else xp.asarray(x0))
                labels.append(dict(arm=arm,history=history,history_sign=1 if history==0 else -1,kicked=bool(kicked)))
    h=xp.stack(states)
    x=xp.stack(positions)
    tangent=xp.asarray(np.tile(np.eye(3),(len(states),P['centers'],1,1)))
    startfields=s.pack(h)
    startx=s.cpu(x)
    energy0=.5*s.inner(h,h)
    budget=xp.zeros(len(states));cflmax=xp.zeros(len(states))
    times,xrows,frows,diagrows,responsefields=[],[],[],[],[]
    t=0.
    for target in np.arange(0,end+1e-9,P['save_dt']):
        while t<target-1e-12:
            step=min(dt,target-t)
            h,x,tangent,db,cfl=s.step(h,x,tangent,step)
            budget+=db;cflmax=xp.maximum(cflmax,cfl);t+=step
        times.append(float(target));xrows.append(s.cpu(x));frows.append(s.cpu(tangent));diagrows.append(s.diagnostics(h,x,tangent,budget,energy0,cflmax))
        # Field responses at declared horizons; differences include all retained modes.
        if any(abs(target-v)<1e-8 for v in P['horizons']+P['long_horizons']):
            responsefields.append((float(target),s.pack(h[1::2]-h[::2])))
    arrays=dict(modes=s.modes_cpu,prep_time=prept,prep_positions=prepx,prep_velocity=prepu,prep_gradient=prepg,
                prep_fields=preph,prep_diagnostics=prepdiag,writer=s.pack(writer[None])[0],probe=s.pack(probe[None])[0],
                start_fields=startfields,start_positions=startx,time=times,positions=xrows,tangent=frows,diagnostics=diagrows,
                end_fields=s.pack(h),response_time=[v[0] for v in responsefields],response_fields=[v[1] for v in responsefields])
    np.savez_compressed(raw_path,**arrays)
    maxima={name:float(np.max(abs(np.asarray(prepdiag)[...,i]))) if i==1 else float(np.max(np.asarray(prepdiag)[...,i])) for i,name in enumerate(DIAG_NAMES)}
    for i,name in enumerate(DIAG_NAMES):
        values=np.asarray(diagrows)[...,i]
        maxima[name]=max(maxima[name],float(np.max(abs(values))) if i==1 else float(np.max(values)))
    template=[]
    for arm,state in variants.items():
        template.append(dict(arm=arm,projection=s.cpu(s.inner(state,writer)).tolist()))
    meta=dict(config=config,status='complete',sources=sources,protocol_sha256=sources['protocol.json'],raw_sha256=digest(raw_path),
              branches=labels,diagnostic_columns=DIAG_NAMES,maxima=maxima,intervention_checks=checks,template_projection=template,
              seconds=time.time()-tic,python=sys.version,platform=platform.platform(),numpy=np.__version__,
              marker_meaning='10 centers and six passive finite-neighbor probes each; fixed geometry for matched causal arms; advected native histories for prediction',
              preparation_count=2,post_branch_count=len(labels),settings_frozen_before_execution=True)
    if backend=='gpu':
        meta.update(cupy=xp.__version__,device=xp.cuda.runtime.getDeviceProperties(0)['name'].decode())
    dump(meta_path,meta)
    print(f'{stem} complete {meta["seconds"]:.1f}s; tail {maxima["tail_enstrophy"]:.3g}; volume {maxima["tangent_volume_error"]:.3g}',flush=True)
    if backend=='gpu':
        xp.get_default_memory_pool().free_all_blocks()
    return meta


def configs(stage):
    base=P['train_seeds']+P['validation_seeds']+P['test_seeds']
    if stage=='pilot':
        return [dict(seed=seed,n=n) for seed in P['pilot_seeds'] for n in P['grids']]
    if stage=='primary':
        return [dict(seed=seed) for seed in base]
    if stage=='refinement':
        return ([dict(seed=seed,n=n) for seed in P['test_seeds'] for n in P['grids'][1:]]+
                [dict(seed=seed,dt=P['dt']/2) for seed in P['test_seeds']]+
                [dict(seed=seed,n=P['grids'][-1],dt=P['dt']/2) for seed in P['test_seeds'][:2]])
    if stage=='sensitivity':
        return ([dict(seed=seed,end=max(P['long_horizons'])) for seed in P['long_horizon_seeds']]+
                [dict(seed=seed,amplitude=P['probe_velocity_rms']/2) for seed in P['amplitude_control_seeds']]+
                [dict(seed=seed,nu=nu) for seed in P['viscosity_control_seeds'] for nu in P['viscosity_controls']])
    raise ValueError(stage)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['pilot','primary','refinement','sensitivity'],default='pilot')
    ap.add_argument('--backend',choices=['cpu','gpu'],default='gpu');ap.add_argument('--out',type=Path);ap.add_argument('--seed',type=int)
    ap.add_argument('--n',type=int);ap.add_argument('--dt',type=float);ap.add_argument('--nu',type=float);ap.add_argument('--end',type=float)
    args=ap.parse_args()
    tasks=configs(args.stage) if args.seed is None else [dict(seed=args.seed,n=args.n,dt=args.dt,nu=args.nu,end=args.end)]
    execution=[]
    for task in tasks:
        execution.append(run_one(**task,backend=args.backend,out=args.out))
        folder=ROOT/'data' if args.out is None else args.out
        dump(folder/('execution-'+args.stage+'-'+args.backend+'.json'),dict(stage=args.stage,completed=len(execution),planned=len(tasks),configs=[m['config'] for m in execution],sources={f:digest(ROOT/f) for f in ['solver.py','run.py','protocol.json']}))
