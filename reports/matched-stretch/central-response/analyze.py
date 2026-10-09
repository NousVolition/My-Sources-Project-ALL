"""Read saved unforced fields; never evolve or alter a simulation.

Run from the original workspace with work/study-env/Scripts/python.exe.
The public measurements.json is sufficient for rebuilding the charts.
"""
import os, sys, json, hashlib, argparse
from pathlib import Path
from datetime import datetime, timezone
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'mpl-cache'))
sys.path.insert(0, str(ROOT.parents[1]/'adversarial-vortex-study'))
sys.path.insert(0, str(ROOT.parent))
from solver import Solver
from numerics import Flow

def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def save(p, d): p.write_text(json.dumps(d, indent=2, allow_nan=False)+'\n', encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(2**20), b''): h.update(chunk)
    return h.hexdigest()
def rel(a,b): return float(abs(a-b)/max(abs(b), 1e-30))

def chord(s, mag, point, direction, height):
    """Half-local-spin chords through a moving marker, not an identified surface."""
    d=direction/np.linalg.norm(direction)
    e=np.cross(d,np.eye(3)[np.argmin(abs(d))]); e/=np.linalg.norm(e)
    f=np.cross(d,e); ds=s.dx/4; distances=np.arange(0,3+ds/2,ds)
    widths=[]; censored=0
    for theta in np.arange(12)*np.pi/12:
        axis=np.cos(theta)*e+np.sin(theta)*f; halves=[]
        for sign in (-1,1):
            points=point+sign*distances[:,None]*axis
            values=s.sample(mag[None],points)[:,0]
            hit=np.flatnonzero(values<height/2)
            if not len(hit): censored+=1; break
            k=int(hit[0])
            if k==0: halves.append(0.); continue
            fraction=(values[k-1]-height/2)/(values[k-1]-values[k])
            halves.append(float(distances[k-1]+fraction*ds))
        if len(halves)==2: widths.append(sum(halves))
    return {'minimum':min(widths) if widths else None,
            'minimum_cells':min(widths)/s.dx if widths else None,
            'censored_directions':censored,
            'definition':'Half of interpolated scalar spin at the selected material marker; 12 transverse directions; not an automatic vortex boundary.'}

def extract():
    ROOT.mkdir(parents=True,exist_ok=True)
    adv=ROOT.parents[1]/'adversarial-vortex-study'
    proto=read(adv/'protocol.json')
    current_hashes={k:sha(adv/k) for k in ('solver.py','initial_design.py','protocol.py','run_study.py')}
    result={'created':datetime.now(timezone.utc).isoformat(),'external_force':0,'viscosity':.001,
            'original':{},'surroundings':{},'local_terms':{},'verification':[],'source_hashes':current_hashes,
            'scope':'Saved-data analysis of two distinct starts. Original matched-stretch uses fixed regions. Separate periodic surroundings study includes material markers. No new evolution.'}
    for p in sorted((adv/'runs').glob('*/result.json')):
        d=read(p)
        if d.get('status')!='complete':continue
        assert d['source_hashes']==current_hashes
        assert d['job']['nu']==.001
        result['surroundings'][p.parent.name]={'job':d['job'],'rows':d['rows'],'source_result_sha256':sha(p)}
    # Freeze the completed-record list: active records never enter this analysis.
    save(ROOT/'input-snapshot.json',{'completed':list(result['surroundings']), 'created':result['created']})
    for suffix in ('base','half'):
        folder=ROOT.parent/'runs'/('baseline-n64-'+suffix)
        p=folder/'result.json';d=read(p); assert d['status']=='complete'
        f=Flow(64,.001,workers=1);a=np.arange(64)*f.dx-3
        x,y,z=np.meshgrid(a,a,a,indexing='ij',sparse=True)
        central=(x*x+y*y<=.4**2)&(abs(z)<=2.5)
        seam=(abs(x)>=2.625)|(abs(y)>=2.625)|(abs(z)>=2.625)
        rows=[]
        for index,r in enumerate(d['series']):
            path=folder/f'field-{index:03d}.npy'; digest=sha(path)
            h=np.load(path,allow_pickle=False);assert np.isfinite(h).all()
            # Independent NumPy inverse FFT for original stored W and energy.
            w=np.fft.irfftn(f.curl(h),s=(64,)*3,axes=(-3,-2,-1));mag=np.linalg.norm(w,axis=0)
            u=np.fft.irfftn(h,s=(64,)*3,axes=(-3,-2,-1));energy=float(.5*np.sum(u*u)*f.dx**3)
            assert rel(float(mag.max()),r['Wmax'])<1e-12 and rel(energy,r['energy'])<1e-12
            ix=np.unravel_index(mag.argmax(),mag.shape)
            cix=np.unravel_index(np.where(central,mag,-1).argmax(),mag.shape)
            row={k:r[k] for k in ('t','Wmax','high_band_enstrophy_fraction','width_at_global_peak')}
            row.update(central_W=float(mag[central].max()),join_band_W=float(mag[seam].max()),
                peak_location=[float(a[i]) for i in ix],central_peak_location=[float(a[i]) for i in cix],
                peak_in_central_region=bool(central[ix]),peak_in_join_band=bool(seam[ix]))
            rows.append(row)
            result['verification'].append({'run':folder.name,'file':path.name,'sha256':digest,'W_energy_match':True})
        result['original'][suffix]={'rows':rows,'dt':d['dt'],'source_result_sha256':sha(p)}
        print('Read original '+suffix+' '+str(len(rows))+' fields',flush=True)
    # Local rates at stored fields every 0.10. Curves above retain 0.02 observations.
    selected=[]
    for case in ('aligned','compressive'):
        for n,suffix in ((112,'base'),(112,'half'),(160,'base')):
            jid=f'{case}-fourier-n{n}-{suffix}'
            if jid in result['surroundings']:selected.append(jid)
    for jid in selected:
        d=result['surroundings'][jid];job=d['job'];s=Solver(job['n'],'fourier',job['nu'],workers=1)
        rows=[];reference={round(r['t'],8):r for r in d['rows']}
        for path in sorted((adv/'runs'/jid).glob('field-t*.npz')):
            digest=sha(path)
            with np.load(path,allow_pickle=False) as z:
                h=z['h'];points=z['markers'][:64];t=float(z['t'])
                assert json.loads(str(z['source_hashes']))==current_hashes
                assert json.loads(str(z['job_json']))==job
            assert np.isfinite(h).all() and np.isfinite(points).all()
            r=reference[round(t,8)];wh=s.curl(h);w=s.real(wh);mag=np.linalg.norm(w,axis=0)
            wm=s.sample(w,points);norm2=np.sum(wm*wm,axis=1);norm=np.sqrt(norm2)
            u=s.real(h);um=s.sample(u,points)
            gu=np.zeros((64,3,3));gw=np.zeros((64,3,3))
            for i in range(3):
                gu[:,i,:]=s.sample(np.stack([s.real(1j*s.k[j]*h[i]) for j in range(3)]),points)
                gw[:,i,:]=s.sample(np.stack([s.real(1j*s.k[j]*wh[i]) for j in range(3)]),points)
            st=np.einsum('nij,nj->ni',gu,wm)
            vi=s.sample(s.real(-s.nu*s.k2*wh),points)
            rh=s.rhs(h)[0]
            dw=s.sample(s.real(s.curl(rh)),points)+np.einsum('nij,nj->ni',gw,um)
            residual=dw-st-vi
            def rate(v): return np.sum(wm*v,axis=1)/np.maximum(norm2,1e-30)
            a,b,c=rate(st),rate(vi),rate(residual)
            weight=norm2/norm2.sum()
            weighted=lambda v:float(np.sum(v*weight))
            assert rel(float(mag.max()),r['Wmax'])<1e-12
            assert rel(float(norm.max()),r['central_marker_Wmax'])<1e-12
            error=abs(weighted(a)-r['central_marker_parallel_stretch_spin_weighted'])
            assert error<1e-10
            best=int(norm.argmax());height=float(s.sample(mag[None],points[best:best+1])[0,0])
            width=chord(s,mag,points[best],wm[best],height)
            peak=np.array(r['peak_location']);delta=(points-peak+3)%6-3
            row={'t':t,'weighted_stretch':weighted(a),'weighted_viscous':weighted(b),
                 'weighted_discrete_interpolation_residual':weighted(c),'weighted_total':weighted(rate(dw)),
                 'marker_spin':norm.tolist(),'marker_locations':points.tolist(),
                 'marker_stretch':a.tolist(),'marker_viscous':b.tolist(),
                 'marker_discrete_interpolation_residual':c.tolist(),
                 'marker_max_index':best,'marker_max_spin':float(norm[best]),
                 'initial_origin_marker_spin':float(norm[32]),
                 'half_local_spin_width':width,
                 'global_peak_distance_to_nearest_central_marker':float(np.linalg.norm(delta,axis=1).min()),
                 'global_peak_location':peak.tolist()}
            assert abs(row['weighted_total']-sum(row[k] for k in ('weighted_stretch','weighted_viscous','weighted_discrete_interpolation_residual')))<1e-10
            rows.append(row)
            result['verification'].append({'run':jid,'file':path.name,'sha256':digest,'W_marker_stretch_match':True})
            print(jid+' t='+str(round(t,2))+' checked',flush=True)
        result['local_terms'][jid]=rows
        save(ROOT/'measurements.partial.json',result)
    result['notes']=[
        'Initial construction, grid, integrator and domain differ between the original matched-stretch and surroundings studies; do not compare their amplitudes as replicas.',
        'Markers label moving fluid points, not vortex identity or the boundary of a material volume. They start on a centreline and do not sample the whole core.',
        'Weighted rates are sum(|omega|^2 * local fractional rate)/sum(|omega|^2). They are not the derivative of a maximum or of the weighted average.',
        'Viscous contribution may be positive at a local point. No constant-brake interpretation is assumed.',
        'The discrete/interpolation residual compares the saved Fourier RHS plus material advection with the continuum stretching-plus-viscosity decomposition sampled at markers.',
        'Physical width is a half-local-spin chord at the currently strongest sampled marker. The marker may change between outputs; threshold changes with local spin.',
        'Exodus data are omitted unless a complete saved run exists; no new run is launched by this analysis.'
    ]
    save(ROOT/'measurements.json',result)
    print(json.dumps({'complete_records':len(result['surroundings']),'verified_fields':len(result['verification']),'new_simulations':0}),flush=True)

if __name__=='__main__': extract()
