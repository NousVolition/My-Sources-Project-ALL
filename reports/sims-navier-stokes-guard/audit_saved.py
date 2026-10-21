"""Read-only numerical audit of saved arrays; never integrates a trajectory."""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import json, hashlib
import numpy as np
from scipy.fft import fft2, ifft2
R=Path(__file__).resolve().parent
read=lambda n:json.loads((R/n).read_text(encoding='utf-8'))
p=read('protocol.json'); cs=read('results/cases.json'); summary=read('results/summary.json')
lookup={(c['mode'],c['weight'],c['pulse']):c['id'] for c in cs}
base=np.array([lookup[c['mode'],c['weight'],'none'] for c in cs])
checks=[]
def check(name,value,limit):checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
for name,expected in read('freeze.json')['files'].items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==expected,name
for label in summary['batch_labels']:
    d=read('results/'+label+'.json');z=np.load(R/'results'/f'{label}.npz');f=z['fields'];ts=z['time'];n=d['n']
    assert f.shape==(5,45,6,n,n) and np.isfinite(f).all()
    assert np.array_equal(ts,np.array(p['time']['field_times']))
    rec=[next(r for r in d['records'] if abs(r['t']-t)<1e-12) for t in ts]
    x,y=np.meshgrid(np.arange(n)*2*np.pi/n,np.arange(n)*2*np.pi/n,indexing='ij')
    expected=np.stack([np.sin(x)*np.cos(y)-.25*np.sin(2*x+y)-.3*np.cos(x-2*y),
        -np.cos(x)*np.sin(y)+.5*np.sin(2*x+y)-.15*np.cos(x-2*y),
        -.55+.15*np.cos(x)+.1*np.sin(y),.35+.1*np.sin(x)-.15*np.cos(y),
        np.full_like(x,.5),np.full_like(x,.5)])
    check(label+'_analytical_start',np.max(np.abs(f[0]-expected)),1e-12)
    actual={
        'energy':.5*np.mean(f[:,:,0]**2+f[:,:,1]**2,axis=(-2,-1)),
        'mean_velocity':np.mean(f[:,:,:2],axis=(-2,-1)),
        'mean_guard':np.mean(f[:,:,4:6],axis=(-2,-1)),
        'gap_rms':np.sqrt(np.mean((f[:,:,3]-f[:,:,2])**2,axis=(-2,-1))),
        'response_rms':np.sqrt(np.mean((f-f[:,base])**2,axis=(-2,-1))),
        'flow_response_rms':np.sqrt(np.mean(np.sum((f[:,:,:2]-f[:,base,:2])**2,axis=2),axis=(-2,-1)))
    }
    k=np.rint(np.fft.fftfreq(n)*n).astype(int);hat=fft2(f,axes=(-2,-1),norm='forward',workers=1)
    dx=ifft2(1j*k[:,None]*hat,axes=(-2,-1),norm='forward',workers=1).real
    dy=ifft2(1j*k[None,:]*hat,axes=(-2,-1),norm='forward',workers=1).real
    check(label+'_saved_divergence',np.max(np.abs(dx[:,:,0]+dy[:,:,1])),1e-10)
    nu=.02+.04*(f[:,:,4]+f[:,:,5])
    for c in cs:
        if c['mode']=='no_feedback':nu[:,c['id']]=.06
    strain=dx[:,:,0]**2+dy[:,:,1]**2+.5*(dy[:,:,0]+dx[:,:,1])**2
    actual['dissipation']=np.mean(2*nu*strain,axis=(-2,-1))
    actual['enstrophy']=.5*np.mean((dx[:,:,1]-dy[:,:,0])**2,axis=(-2,-1))
    actual['mean_viscosity']=np.mean(nu,axis=(-2,-1))
    for key,value in actual.items():check(label+'_saved_'+key,np.max(np.abs(value-np.array([r[key] for r in rec]))),1e-12)
    energy=np.array([r['energy'] for r in d['records']]);loss=np.array([r['loss'] for r in d['records']])
    check(label+'_energy_does_not_increase',max(0,np.diff(energy,axis=0).max()),1e-12)
    check(label+'_integrated_dissipation_nondecreasing',max(0,-np.diff(loss,axis=0).min()),1e-12)
    check(label+'_saved_budget_identity',np.max(np.abs((energy-energy[0]+loss)/energy[0])),p['hard_limits']['energy_relative'])
    check(label+'_batch_time',d['elapsed_seconds'],p['resources']['max_batch_seconds'])
check('total_integration_time',summary['elapsed_seconds'],p['resources']['max_total_seconds'])
check('result_storage_GiB',sum(f.stat().st_size for f in (R/'results').iterdir() if f.is_file())/2**30,p['resources']['max_storage_GiB'])
out=dict(passed=all(c['passed'] for c in checks),frozen_hashes_matched=True,checks=checks,
         scope='Independent remeasurement of saved fields and budgets, not additional integration. The failed spatial refinement screens in verification.json remain in force.')
(R/'results/saved-field-audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(dict(passed=out['passed'],checks=len(checks),failures=[c for c in checks if not c['passed']]),indent=2))
