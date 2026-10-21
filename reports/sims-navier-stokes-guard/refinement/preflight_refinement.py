"""Check exact source reuse, common starting fields and feasible step sizes."""
from pathlib import Path
import os,json,hashlib
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np
from scipy.fft import fft2,ifft2
from threadpoolctl import threadpool_info
from model import Model,cases
R=Path(__file__).resolve().parent
p=json.loads((R/'protocol.json').read_text());B=Path(p['refinement']['parent_directory']);cs=cases(p)
assert cs==json.loads((B/'results/cases.json').read_text())
assert hashlib.sha256((R/'model.py').read_bytes()).hexdigest()==p['refinement']['parent_model_sha256']
checks=[]
def check(name,value,limit):
    row=dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit));checks.append(row)
    assert row['passed'],row
m=Model(64,p,cs);h=m.start();f=m.real(h)
old=np.load(B/'results/n48-dt0.005.npz')['fields'][0]
k=np.rint(np.fft.fftfreq(48)*48).astype(int)%64
hat=np.zeros((45,6,64,64),dtype=complex);hat[...,k[:,None],k[None,:]]=fft2(old,axes=(-2,-1),norm='forward')
check('same_finite_start_48_to_64',np.max(abs(f-ifft2(hat,axes=(-2,-1),norm='forward').real)),1e-12)
q=Model(32,p,cs)
check('same_pulse_profile_32_to_64',np.max(abs(q.profile-m.profile[::2,::2])),1e-12)
rhs,aux=m.rhs(0,h);df=m.real(rhs);d=m.diagnostics(h)
check('64_grid_discrete_energy_identity',np.max(abs(np.mean(f[:,0]*df[:,0]+f[:,1]*df[:,1],axis=(-2,-1))+d['dissipation'])),1e-12)
check('64_grid_initial_divergence',aux['divergence'],p['hard_limits']['divergence'])
check('initial_CFL_larger_step',aux['speed']*.005*64/(2*np.pi),p['hard_limits']['CFL'])
check('diffusive_RK4_scale_estimate',.1*2*21**2*.005,2.785)
pool=threadpool_info();assert all(x['num_threads']==1 for x in pool),pool
out=dict(passed=True,checks=checks,blas_pools=pool,fft_workers=1,initial_fourier_bytes_sha256=hashlib.sha256(h.tobytes()).hexdigest(),scope='Operator, starting-field and resource checks only; full trajectories follow.')
(R/'preflight.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(dict(passed=True,checks=checks,threads=[x['num_threads'] for x in pool]),indent=2))
