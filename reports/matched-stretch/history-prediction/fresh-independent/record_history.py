"""Full parent joint-Heun histories, with exclusive ownership and saved recovery state."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import sys
import json
import time
import argparse
import datetime
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(BASE))
from simulate import Flow, L, initial, sha, save_json
from sampler import sample
sys.path.insert(0, str(BASE/'resolution-followup'))
from analyze import source_hash_matches


def sources():
    paths = [HERE/'record_history.py', HERE/'protocol.json', HERE/'sampler.py', HERE/'fresh_features.py',
             BASE/'simulate.py', BASE.parent/'numerics.py', BASE/'features.py', BASE/'evaluate.py', BASE/'protocol.json']
    return {str(p.relative_to(BASE.parent)).replace('\\', '/'): sha(p) for p in paths}


def name_for(seed, config):
    return f"s{seed}_n{config['n']}_nu{config['nu']}_dt{config['dt']}"


def validate(path, config=None):
    meta = json.loads(path.with_suffix('.json').read_text())
    if config is not None and meta['config'] != config:
        raise ValueError('Configuration mismatch '+str(path))
    if sha(path) != meta['sha256']:
        raise ValueError('Archive checksum mismatch '+str(path))
    for rel, digest in meta['sources'].items():
        p = BASE.parent/'numerics.py' if rel == 'numerics.py' else BASE/rel
        if not source_hash_matches(p.read_bytes(), digest):
            raise ValueError('Changed parent source '+rel)
    if 'study_sources' in meta:
        if set(meta['study_sources'])!=set(sources()):raise ValueError('Changed executed study source set')
        for rel,digest in meta['study_sources'].items():
            if not source_hash_matches((BASE.parent/rel).read_bytes(),digest):raise ValueError('Changed executed source '+rel)
    c = meta['config']
    if meta['max_cfl'] >= .5 or max(abs(d['energy_budget_residual']) for d in meta['diagnostics']) >= .005:
        raise ValueError('CFL/energy check failed')
    if max(d['divergence_max'] for d in meta['diagnostics']) >= 1e-10:
        raise ValueError('Divergence check failed')
    with np.load(path, allow_pickle=False) as z:
        if not all(np.isfinite(z[k]).all() for k in z.files):
            raise ValueError('Nonfinite archive')
        np.testing.assert_allclose(z['time'], np.arange(round(c['end']/c['save_dt'])+1)*c['save_dt'], rtol=0, atol=3e-16)
        expected = np.random.default_rng(meta['seed']+80000).uniform(-L/2,L/2,(c['markers'],3))
        np.testing.assert_array_equal(z['positions'][0], expected)
        np.testing.assert_array_equal(z['tangent'][0], np.broadcast_to(np.eye(3),(c['markers'],3,3)))
    return meta


def run(seed, config, destination=None):
    p=json.loads((HERE/'protocol.json').read_text())
    if seed in p['test_seeds']:
        from fit_models import load_models
        lock=json.loads((HERE/'test-lock.json').read_text())
        manifest=HERE/lock['manifest']
        if sha(manifest)!=lock['manifest_sha256']:raise ValueError('Test model lock changed')
        load_models(manifest)
        publication=json.loads((HERE/'model-lock-publication.json').read_text())
        if not publication['checks_success'] or publication['test_lock_sha256']!=sha(HERE/'test-lock.json'):
            raise ValueError('Test generation requires verified model-lock publication')
    started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()
    destination = Path(destination) if destination else HERE/'recorded-data'
    destination.mkdir(parents=True, exist_ok=True)
    name = name_for(seed, config); path = destination/(name+'.npz')
    if path.exists():
        validate(path, config)
        return path
    claims = destination.parent/'claims'; claims.mkdir(exist_ok=True)
    claim = claims/(name+'.json'); frozen = sources()
    with claim.open('x') as stream:
        json.dump({'name': name, 'pid': os.getpid(), 'status': 'running', 'sources': frozen}, stream)
    recovery = destination.parent/'recovery'; recovery.mkdir(exist_ok=True)
    checkpoint = recovery/(name+'_latest.npz')
    if checkpoint.exists():
        raise ValueError('Existing partial checkpoint preserved; inspect before restart')
    f = Flow(config['n'], config['nu'], workers=1); h = initial(f, seed)
    x = np.random.default_rng(seed+80000).uniform(-L/2,L/2,(config['markers'],3))
    M = np.broadcast_to(np.eye(3),(len(x),3,3)).copy()
    dt=config['dt']; nstep=round(config['end']/dt); stride=round(config['save_dt']/dt)
    assert np.isclose(nstep*dt,config['end']) and np.isclose(stride*dt,config['save_dt'])
    arrays={k:[] for k in ['time','positions','velocity','gradient','tangent']}
    diag=[]; start=time.perf_counter(); energy0=f.inner(h,h)/2
    max_cfl=0.; integrated_diss=0.
    try:
        for step in range(nstep+1):
            u=f.real(h); vel,J=sample(u,x)
            if step % stride == 0:
                for key,value in zip(arrays,[step*dt,x.copy(),vel,J,M.copy()]): arrays[key].append(value)
                wh=f.curl(h)
                divergence=f.real(1j*sum(k*c for k,c in zip(f.k,h)))
                energy=f.inner(h,h)/2
                diag.append({'t':step*dt,'energy':energy,'energy_budget_residual':(energy-energy0+integrated_diss)/energy0,
                    'divergence_max':float(abs(divergence).max()),
                    'high_band_energy_fraction':f.inner(h*f.high,h*f.high)/max(2*energy,1e-30),
                    'max_tangent_condition':float(np.linalg.cond(M).max()),'max_det_error':float(abs(np.linalg.det(M)-1).max())})
                if step % round(.1/dt) == 0 or step == nstep:
                    if sources() != frozen: raise ValueError('Active numerical/feature/model source changed')
                    temp=checkpoint.with_suffix('.tmp.npz')
                    np.savez_compressed(temp,field=h,positions=x,tangent=M,step=step,
                        integrated_diss=integrated_diss,max_cfl=max_cfl,
                        **{'saved_'+key:np.asarray(value) for key,value in arrays.items()})
                    temp.replace(checkpoint)
                    save_json(checkpoint.with_suffix('.json'), {'seed':seed,'config':config,'sources':frozen,
                        'checkpoint_sha256':sha(checkpoint),'diagnostics':diag,'elapsed':step*dt})
                    print(name+f' saved t={step*dt:.3f}',flush=True)
            if step == nstep: break
            # Exactly the original parent joint Heun arithmetic and ordering.
            r1,a,speed1=f.rhs(h); hp=h+dt*r1; xp=x+dt*vel
            dM=J@M; Mp=M+dt*dM
            vp,Jp=sample(f.real(hp),xp); r2,b,speed2=f.rhs(hp)
            h=f.project(h+dt/2*(r1+r2)); x=x+dt/2*(vel+vp); M=M+dt/2*(dM+Jp@Mp)
            integrated_diss+=dt/2*(a[1]+b[1]); max_cfl=max(max_cfl,max(speed1,speed2)*dt/f.dx)
            if max_cfl > .5 or not np.isfinite(M).all() or not np.isfinite(h).all():
                raise FloatingPointError('Parent numerical safety check failed')
        np.savez_compressed(path,**arrays)
        meta={'name':name,'independent_group':seed,'seed':seed,'config':config,
              'sources':{p.name:sha(p) for p in [BASE/'simulate.py',BASE.parent/'numerics.py',BASE/'protocol.json']},
              'study_sources':frozen,'sha256':sha(path),'seconds':time.perf_counter()-start,
              'max_cfl':max_cfl,'diagnostics':diag,'Re_Urms_L_over_nu':L/config['nu'],
              'started_UTC':started_UTC,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
        save_json(path.with_suffix('.json'),meta); validate(path,config)
        checkpoint.unlink(); checkpoint.with_suffix('.json').unlink()
        save_json(claim,{'name':name,'pid':os.getpid(),'status':'complete','sources':frozen,'sha256':meta['sha256']})
        print(name+' complete',flush=True)
        return path
    except Exception as exc:
        save_json(claim,{'name':name,'pid':os.getpid(),'status':'failed','sources':frozen,'error':repr(exc)})
        raise


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--grid',type=int,required=True); ap.add_argument('--dt',type=float); ap.add_argument('--nu',type=float)
    args=ap.parse_args(); p=json.loads((HERE/'protocol.json').read_text())
    run(args.seed,{**p['base'],'n':args.grid,'dt':args.dt or p['base']['dt'],'nu':args.nu or p['base']['nu']})


if __name__ == '__main__': main()
