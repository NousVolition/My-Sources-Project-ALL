"""Separate initial-field diagnostic; no time evolution or changes to running solvers."""
import gc,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source'))
from numerics import Flow,L
def save(name,value):(ROOT/name).write_bytes((json.dumps(value,indent=2,allow_nan=False)+'\n').encode())
def tube(f):
    modes=np.rint(np.fft.fftfreq(f.n)*f.n).astype(int)
    mx,my=modes[:,None],modes[None,:];k=2*np.pi/L
    ph=np.zeros((f.n,f.n,f.n//2+1),complex)
    ph[:,:,0]=f.n**3*.4*np.pi*.2**2/L**2*np.exp(-.2**2*k*k*(mx*mx+my*my)/4)*(-1.)**(mx+my)
    ph*=f.keep
    return np.stack((1j*f.k[1]*ph,-1j*f.k[0]*ph,np.zeros_like(ph)))
def background(f):
    k=2*np.pi/L;a=np.arange(f.n)*f.dx-L/2
    xyz=np.meshgrid(a,a,a,indexing='ij',sparse=True)
    coords=[np.sin(k*x)/k for x in xyz]
    r2=sum(2*(1-np.cos(k*x))/k**2 for x in xyz)
    env=np.exp(-.04*r2)
    raw=np.stack([c*x*env for c,x in zip((-40.,-40.,80.),coords)])
    return f.project(f.hat(raw))
def measure(f,h):
    wh=f.curl(h);w=f.real(wh);mag=np.sqrt(np.sum(w*w,axis=0))
    idx=np.unravel_index(int(mag.argmax()),mag.shape);v=w[(slice(None),)+idx]
    axis=np.arange(f.n)*f.dx-L/2
    div=float(np.max(abs(f.real(1j*sum(k*x for k,x in zip(f.k,h))))))
    gradient=np.array([[float(f.real(1j*f.k[j]*h[i])[f.n//2,f.n//2,f.n//2]) for j in range(3)] for i in range(3)])
    return dict(W=float(mag[idx]),peak_index=list(map(int,idx)),peak_location=[float(axis[i]) for i in idx],
      width=f.width(mag,idx,v),energy=.5*f.inner(h,h),enstrophy=.5*f.inner(wh,wh),
      high_band_energy=f.inner(h*f.high,h*f.high)/f.inner(h,h),
      high_band_enstrophy=f.inner(wh*f.high,wh*f.high)/f.inner(wh,wh),
      mean=(h[:,0,0,0].real/f.n**3).tolist(),divergence_max=div,center_velocity_gradient=gradient.tolist())
def gate(r):
    reasons=[]
    if (r['width']['minimum_chord_cells'] or 0)<6:reasons.append('width below 6 cells')
    if r['high_band_energy']>.001:reasons.append('high-band energy above 0.1%')
    if r['high_band_enstrophy']>.01:reasons.append('high-band enstrophy above 1%')
    return dict(passed=not reasons,reasons=reasons)
def main():
    reference=json.loads((ROOT/'reference-mean.json').read_text())
    mean=np.array(reference['mean_velocity'])
    protocol=dict(grids=[64,80,112,128,256],L=L,nu=.001,mode='initialization only',
      original='Unchanged supplied start, measured separately; no old checkpoints replaced.',
      candidate='Periodic sine coordinates sin(k*x)/k and periodic squared distance 2*(1-cos(k*x))/k^2 replace raw strain coordinates and envelope; k=2*pi/L. Exact periodized-Gaussian tube coefficients.',
      reference_mean=mean.tolist(),mean_rule='Same fixed mean on every candidate grid, taken from original 256 initialization.',
      energy_rule='No energy normalization; record energy change explicitly.',
      amplitude_rule='Raw center linear strain coefficients remain -40,-40,80; projection can change actual center strain, which is measured.',
      source_sha256=hashlib.sha256((ROOT/'source/numerics.py').read_bytes()).hexdigest())
    save('protocol.json',protocol);rows=[];slices=[]
    for n in protocol['grids']:
        f=Flow(n,workers=2);old=f.initial();newtube=tube(f);newbg=background(f);new=newtube+newbg;new[:,0,0,0]=mean*n**3
        om=measure(f,old);nm=measure(f,new);om['gate']=gate(om);nm['gate']=gate(nm)
        oldbg=f.initial(spin_factor=0.);oldtube=old-oldbg
        oldidx=tuple(om['peak_index']);newidx=tuple(nm['peak_index'])
        contributions={}
        for label,h,idx in [('original_background',oldbg,oldidx),('original_tube',oldtube,oldidx),('candidate_background',newbg,newidx),('candidate_tube',newtube,newidx)]:
            w=f.real(f.curl(h));contributions[label]=w[(slice(None),)+idx].tolist()
        a=np.array(contributions['original_background'])+np.array(contributions['original_tube'])
        b=np.array(contributions['candidate_background'])+np.array(contributions['candidate_tube'])
        assert abs(np.linalg.norm(a)-om['W'])<1e-9 and abs(np.linalg.norm(b)-nm['W'])<1e-9
        error=np.sqrt(f.inner(newtube-oldtube,newtube-oldtube)/f.inner(oldtube,oldtube))
        row=dict(n=n,original=om,candidate=nm,energy_ratio=nm['energy']/om['energy'],tube_relative_L2_difference=float(error),contributions_at_each_combined_peak=contributions)
        rows.append(row)
        if n in (128,256):
            for label,h in [('original',old),('candidate',new)]:
                w=f.real(f.curl(h));mag=np.sqrt(np.sum(w*w,axis=0))
                slices.append(dict(n=n,kind=label,z=0.,mag=mag[:,:,n//2].tolist()))
        save('results.json',dict(status='running' if n!=256 else 'complete',new_time_evolution=False,rows=rows));save('slices.json',slices)
        print(json.dumps(dict(n=n,original_width=om['width']['minimum_chord_cells'],candidate_width=nm['width']['minimum_chord_cells'],candidate_W=nm['W'],candidate_gate=nm['gate'],energy_ratio=row['energy_ratio'])),flush=True)
        del f,old,newtube,newbg,new,oldbg,oldtube,w;gc.collect()
if __name__=='__main__':main()
