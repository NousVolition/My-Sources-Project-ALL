import numpy as np

NAMES=['energy','dissipation','helicity','speed_fourth_moment','flux_k4','energy_fraction_forced','energy_fraction_k_ge_4','integral_length',
       'local_u_x','local_u_y','local_u_z']+[f'local_gradient_{i}{j}' for i in 'xyz' for j in 'xyz']+['local_vortex_stretching','local_strain_norm','local_du_dt_x','local_du_dt_y','local_du_dt_z']

def inputs(z):
 """Only recorded times at/before the impulse. Future responses are never read."""
 h=np.asarray(z['history']);a=np.asarray(z['strengths']);events=np.asarray(z['events']);r=np.asarray(z['initial_rates'])
 if h.ndim!=3 or h.shape[2]!=len(NAMES):raise ValueError('Unexpected history shape')
 cur=h[:,0,:];changes=(h[:,1:,:]-cur[:,None,:]).reshape(len(cur),-1)
 c=np.repeat(cur,len(a),axis=0);past=np.repeat(changes,len(a),axis=0)
 intervention=np.c_[np.tile(a,len(events)),np.repeat(events,len(a))]
 conventional=np.c_[c[:,:8],intervention]
 present=np.c_[c,intervention,r.reshape(-1,3)]
 return conventional,present,past

def target(z):return np.trapezoid(z['energy_ratio'],z['response_time'][...,None],axis=1).reshape(-1)

def derange(n,seed):
 rng=np.random.default_rng(seed)
 for _ in range(10000):
  p=rng.permutation(n)
  if np.all(p!=np.arange(n)):return p
 raise RuntimeError('No derangement')

def shuffle(h,n,seed):
 a=np.asarray(h).reshape(n,-1,h.shape[-1]);return a[derange(n,seed)].reshape(h.shape)
