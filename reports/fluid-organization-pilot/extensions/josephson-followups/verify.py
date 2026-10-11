"""Recompute checks from saved arrays and verify frozen models and file provenance."""
import json,hashlib,sys,platform
import numpy as np
import scipy,matplotlib
from support import ROOT,P,save
from prediction import predict,shuffled_history
def main():
    checks=[]
    def check(name,ok,value=None):
        checks.append({'check':name,'pass':bool(ok),'observed':value})
    def read(name):return json.loads((ROOT/f'results/{name}.json').read_text())
    for name in ('long_cycles','capacitance','fluid_prediction'):
        for k,v in read(name)['gates'].items():check(name+'/'+k,v)
    for name,keys in [('ramps',['step_gate_2e_6']),('noise',['resolution_gate','covariance_gate_0_02'])]:
        for k in keys:check(name+'/'+k,read(name)[k])
    l=np.load(ROOT/'data/long_cycles.npz');check('324 cycle checkpoints',l['state'].shape==(324,2))
    check('direct time continuation recomputed',np.max(abs(l['direct_voltage']-l['voltage'][:,l['direct_indices']]))<2e-6)
    check('cycle-energy recomputed',np.max(abs(l['cycle_balance']))<1e-6)
    check('cycle time strictly increases',np.all(np.diff(l['cycle_time'],axis=1)>0))
    n=np.load(ROOT/'data/noise.npz');ix=int(np.argmin(abs(n['time']-100)))
    for name in ('fine','coarse'):
        computed=(n[name+'_phase'][-1]-n[name+'_phase'][ix]).sum(-1)/300
        check(name+' noisy voltage reproduced',np.allclose(computed,n[name+'_voltage'],atol=1e-14))
    check('noise clones identical',np.array_equal(n['fine_phase'][0],n['coarse_phase'][0]))
    # Budget diagnostics are reported, not silently promoted to prespecified gates.
    diagnostics={'noise_fine_budget_max':float(abs(n['potential_budget_residual']).max()),'noise_coarse_budget_max':float(abs(n['coarse_potential_budget_residual']).max())}
    c=np.load(ROOT/'data/capacitance.npz')
    for name in ('base','half_step','double_settle'):
        q=c[name+'_states'];check(name+' capacitor continuation',np.allclose(q[1:,0,...,:2],q[:-1,-1,...,:2],rtol=0,atol=1e-12))
        check(name+' capacitor voltage recomputed',np.allclose(2*(q[:,-1,...,0]-q[:,1,...,0])/100,c[name+'_voltage'],rtol=0,atol=1e-12))
    d=np.load(ROOT/'data/circuit_prediction_inputs.npz');g=d['groups'];test=g>=32
    check('circuit run split counts',[len(np.unique(g[g<24])),len(np.unique(g[(g>=24)&(g<32)])),len(np.unique(g[test]))]==[24,8,16])
    check('circuit outcome recomputed',np.allclose(d['y'],2*(d['end_changed'][:,0]-d['end_unchanged'][:,0])/20,atol=1e-14))
    check('full-state independent solver comparator',read('circuit_prediction')['full_state_physics_max_error']<2e-6)
    fluidtrain=np.load(ROOT/'data/fluid_training.npz');fluidtest=np.load(ROOT/'data/fluid_test_features.npz')
    check('fluid train/val/test independence',not(set(fluidtrain['train_groups'])&set(fluidtest['groups']) or set(fluidtrain['validation_groups'])&set(fluidtest['groups']) or set(fluidtrain['train_groups'])&set(fluidtrain['validation_groups'])))
    # Recreate the fixed split-wise null RNG sequence, including train and validation.
    for name in ('circuit','fluid'):
        rng=np.random.default_rng(43120);transformed={}
        for split in ('train','validation','test'):
            if name=='circuit':
                mask=g<24 if split=='train' else ((g>=24)&(g<32) if split=='validation' else test)
                part={k:d[k][mask] for k in ('current','simple','history','groups','stratum')}
            else:
                part={k:(fluidtest[k] if split=='test' else fluidtrain[split+'_'+k]) for k in ('current','history','groups')}
            transformed[split]={'simple_present':part.get('simple',part['current']),'present':part['current'],'history':np.c_[part['current'],part['history']],'shuffled_history':np.c_[part['current'],shuffled_history(part,rng)]}
        out=np.load(ROOT/f'data/{name}_predictions.npz')
        for model in ('simple_present','present','history','shuffled_history'):
            frozen=np.load(ROOT/f'data/{name}_model_{model}.npz');forecast=predict(frozen,transformed['test'][model])
            check(name+'/'+model+' frozen predictions reproduce',np.allclose(forecast,out[model],rtol=0,atol=1e-8))
            check(name+'/'+model+' train-only scaling',np.allclose(frozen['mean'],transformed['train'][model].mean(0),rtol=0,atol=1e-12))
    for p in sorted((ROOT/'fluid/data/fluid').glob('*.json')):
        meta=json.loads(p.read_text());raw=p.with_suffix('.npz')
        check(p.stem+' raw SHA256',hashlib.sha256(raw.read_bytes()).hexdigest()==meta['sha256'])
        check(p.stem+' exact source provenance',all(hashlib.sha256((ROOT/'fluid'/name).read_bytes()).hexdigest()==sha for name,sha in meta['sources'].items()))
        check(p.stem+' energy/clone/CFL',meta['energy_budget_residual']<1e-5 and meta['burn_budget_residual']<1e-5 and meta['clone_error']==0 and meta['max_cfl']<.5 and meta['divergence_max']<1e-8)
    for p in sorted((ROOT/'data').glob('*.npz')):
        q=np.load(p);check(p.name+' finite arrays',all(np.isfinite(q[k]).all() for k in q.files if np.issubdtype(q[k].dtype,np.number)))
    confirmation=read('fluid_confirmation');z=np.load(ROOT/'data/fluid_confirmation.npz');groups=z['groups']
    check('120 completely new confirmation fields',np.array_equal(np.unique(groups),np.arange(5000,5120)) and np.all(np.unique(groups,return_counts=True)[1]==2))
    check('confirmation excludes old training, validation and test seeds',not(set(groups)&(set(range(1000,1036))|set(range(4000,4012)))))
    locked=json.loads((ROOT/'data/confirmation_frozen_hashes.json').read_text())['files']
    check('confirmation models and protocol remain frozen',all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha for p,sha in locked.items()))
    current=z['current'];history=z['history_features'];shuffled=shuffled_history({'current':current,'history':history,'groups':groups},np.random.default_rng(43132))
    for name,x in [('present',current),('history',np.c_[current,history]),('shuffled_history',np.c_[current,shuffled])]:
        forecast=predict(np.load(ROOT/f'data/fluid_model_{name}.npz'),x)
        check('confirmation '+name+' frozen forecasts reproduce',np.allclose(forecast,z[name],atol=1e-8,rtol=0))
    errors={k:np.sqrt(np.mean((z[k]-z['truth']).reshape(-1,2)**2,axis=1)) for k in ('present','history')};delta=errors['history']-errors['present']
    ix=np.random.default_rng(43131).integers(0,120,(20000,120));ci=np.quantile(delta[ix].mean(1),[.025,.975])
    check('confirmation paired whole-field interval reproduced',np.allclose(ci,confirmation['history_minus_present']['ci95'],atol=1e-12))
    check('confirmation decision follows declared rule',bool(ci[1]<0)==confirmation['primary_improvement_resolved'])
    for k,v in confirmation['gates'].items():check('confirmation/'+k,v)
    save('results/verification.json',{'checks':checks,'passed':sum(c['pass'] for c in checks),'total':len(checks),'all_pass':all(c['pass'] for c in checks),'diagnostics':diagnostics,'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform()}})
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or '.pytest_cache' in p.parts or p.name=='manifest.json':continue
        b=p.read_bytes();files.append({'path':p.relative_to(ROOT).as_posix(),'size':len(b),'sha256':hashlib.sha256(b).hexdigest()})
    save('manifest.json',{'description':'Final source, data and result snapshot; manifest excludes itself. Numerical reproduction may change timing metadata and floating-point last bits.','files':files})
    print(json.dumps({'passed':sum(c['pass'] for c in checks),'total':len(checks),'diagnostics':diagnostics,'files':len(files)}))
    if not all(c['pass'] for c in checks):raise AssertionError([c for c in checks if not c['pass']])
if __name__=='__main__':main()
