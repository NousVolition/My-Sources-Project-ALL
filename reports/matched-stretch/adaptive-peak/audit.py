"""Resolution gate and peak diagnostics from immutable saved original fields."""
import argparse,hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
from scipy import ndimage
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source'))
from numerics import Flow,L

def save(p,v):p.write_bytes((json.dumps(v,indent=2,allow_nan=False)+'\n').encode())
def profile(mag,idx,w,dx):
    n=mag.shape[0];d=w/np.linalg.norm(w);e=np.cross(d,np.eye(3)[np.argmin(abs(d))]);e/=np.linalg.norm(e);f=np.cross(d,e)
    origin=np.array(idx)*dx;radii=np.arange(0,L/2+dx/8,dx/4);level=mag[idx]/2;choices=[]
    for theta in np.linspace(0,np.pi,12,endpoint=False):
        axis=np.cos(theta)*e+np.sin(theta)*f;cross=[]
        for sign in (-1,1):
            points=((origin+sign*radii[:,None]*axis)%L/dx).T
            v=ndimage.map_coordinates(mag,points,order=1,mode='grid-wrap');hits=np.flatnonzero(v<level)
            if len(hits):
                j=int(hits[0]);cross.append(float(radii[j-1]+(level-v[j-1])*(radii[j]-radii[j-1])/(v[j]-v[j-1])))
        if len(cross)==2:choices.append((sum(cross),axis,cross))
    if not choices:return {'censored':True}
    width,axis,cross=min(choices,key=lambda a:a[0]);xi=np.linspace(-3,3,241);s=xi*width/2
    points=((origin+s[:,None]*axis)%L/dx).T
    values=ndimage.map_coordinates(mag,points,order=1,mode='grid-wrap')/mag[idx]
    return dict(censored=False,width=width,width_cells=width/dx,axis=axis.tolist(),half_crossings=cross,scaled_distance=xi.tolist(),scaled_magnitude=values.tolist())

def detailed(f,h,t):
    row=f.observe(h,t);wh=f.curl(h);w=f.real(wh);u=f.real(h);mag=np.sqrt(np.sum(w*w,axis=0));idx=np.unravel_index(int(mag.argmax()),mag.shape);ix=(slice(None),)+idx
    gu=np.array([[f.real(1j*f.k[j]*h[i])[idx] for j in range(3)] for i in range(3)])
    S=(gu+gu.T)*.5;vals,vecs=np.linalg.eigh(S);direction=w[ix]/mag[idx]
    degenerate=vals[-1]-vals[-2]<1e-8*max(1.,abs(vals[-1]))
    angle=None if degenerate else float(np.degrees(np.arccos(np.clip(abs(np.dot(direction,vecs[:,-1])),0,1))))
    transport=np.zeros_like(w);stretch=np.zeros_like(w)
    for i in range(3):
        for j in range(3):
            transport[i]+=u[j]*f.real(1j*f.k[j]*wh[i]);stretch[i]+=w[j]*f.real(1j*f.k[j]*h[i])
    transport=f.real(f.hat(transport)*f.keep);stretch=f.real(f.hat(stretch)*f.keep);visc=f.real(-f.nu*f.k2*wh)
    rh,_,_=f.rhs(h);actual=f.real(f.curl(rh));summed=-transport+stretch+visc
    row['peak_budget']={name:{'vector':v[ix].tolist(),'along_vorticity':float(np.dot(v[ix],direction))} for name,v in [('transport_LHS',transport),('stretching_RHS',stretch),('viscosity_RHS',visc)]}
    row['peak_budget']['actual_time_derivative']=float(np.dot(actual[ix],direction))
    row['peak_budget']['closure_relative_l2']=float(np.linalg.norm(actual-summed)/max(np.linalg.norm(actual),1e-300))
    row['strain_eigenvalues']=vals.tolist();row['alignment_degrees']=angle;row['strain_axis_degenerate']=bool(degenerate)
    row['near_tied_maxima']=int(np.sum(mag>=row['Wmax']*(1-1e-8)))
    row['profile']=profile(mag,idx,w[ix],f.dx)
    radius=np.sqrt(f.k2)/(2*np.pi/L);shell=np.rint(radius).astype(int)
    modalE=.5*f.scale*np.sum(abs(h)**2,axis=0)*f.weights
    modalZ=.5*f.scale*np.sum(abs(wh)**2,axis=0)*f.weights
    counts=np.bincount(shell.ravel(),weights=np.broadcast_to(f.weights,shell.shape).ravel())
    es=np.bincount(shell.ravel(),weights=modalE.ravel());zs=np.bincount(shell.ravel(),weights=modalZ.ravel())
    row['spectra']={'shell_index':list(range(len(es))),'energy_sum':es.tolist(),'enstrophy_sum':zs.tolist(),
                    'energy_per_mode':(es/np.maximum(counts,1)).tolist(),'enstrophy_per_mode':(zs/np.maximum(counts,1)).tolist(),
                    'mode_counts':counts.tolist(),'nominal_component_cutoff':f.n//3,'actual_component_limit':max(float(np.max(abs(k*f.keep)))/(2*np.pi/L) for k in f.k),
                    'actual_radial_limit':float(radius[f.keep].max())}
    return row

def gate(row):
    w=row['width_at_global_peak']['minimum_chord_cells'];reasons=[]
    if w is None or w<6:reasons.append('core narrower than 6 cells')
    if row['high_band_enstrophy_fraction']>.01:reasons.append('high-band enstrophy above 1%')
    if row['high_band_energy_fraction']>.001:reasons.append('high-band energy above 0.1%')
    return {'pass':not reasons,'reasons':reasons}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--study-root',type=Path,default=ROOT.parent);ap.add_argument('--legacy-root',type=Path,default=ROOT.parents[1]/'fixed-cutoff-run');ap.add_argument('--max-128-time',type=float,default=None);args=ap.parse_args()
    assert hashlib.sha256((ROOT/'source/numerics.py').read_bytes()).hexdigest()=='22b169971aaf2efd56d47f685495a575e1622edd56c199b07a785fc97e952f46'
    initial=json.loads((args.study_root/'initial-audit.json').read_text());starts={r['n']:r for r in initial['rows']}
    snapshots=[];runs={};manifest=[]
    for n in (64,128):
        folder=args.study_root/'runs'/f'baseline-n{n}-base';meta=json.loads((folder/'result.json').read_text());series=list(meta['series']);f=Flow(n,workers=2);rows=[]
        if n==128 and args.max_128_time is not None:series=[r for r in series if r['t']<=args.max_128_time+1e-12]
        for i,old in enumerate(series):
            path=folder/f'field-{i:03d}.npy';h=np.load(path);row=detailed(f,h,old['t']);assert abs(row['Wmax']-old['Wmax'])<1e-8
            row['I']=old['I'];row['gate']=gate(row);rows.append(row)
            if n==64 and i in (0,12,24,40):
                w=f.real(f.curl(h));mag=np.sqrt(np.sum(w*w,axis=0));idx=np.unravel_index(mag.argmax(),mag.shape)
                snapshots.append(dict(n=n,t=old['t'],z=-L/2+idx[2]*f.dx,mag=mag[:,:,idx[2]].tolist()))
            manifest.append(dict(run=f'n{n}',time=old['t'],file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        runs[str(n)]=dict(n=n,nu=.001,L=6,dt=meta['dt'],source='Existing fixed matrix; snapshots frozen at audit start',rows=rows)
        print('CHECKED saved grid',n,len(rows),'through',rows[-1]['t'],flush=True)
    n=80;folder=args.legacy_root/'n80';meta=json.loads((folder/'result.json').read_text());f=Flow(n,workers=2);rows=[]
    for old in meta['series']:
        path=folder/old['field'];z=np.load(path);u=np.array([z['u0'],z['u1'],z['u2']]);h=f.hat(u)
        row=detailed(f,h,old['exact_time']);assert abs(row['Wmax']-old['biggest'])<1e-8
        row['I']=old['integral'];row['gate']=gate(row);rows.append(row)
        manifest.append(dict(run='n80',time=old['exact_time'],file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    starts[80]=dict(rows[0],n=80);runs['80']=dict(n=80,nu=.001,L=6,dt=meta['dt'],source='Completed original replay',rows=rows)
    f=Flow(112,workers=2);h=f.initial();r=detailed(f,h,0.);starts[112]=dict(r,n=112)
    print('CHECKED new 112 initialization only',r['width_at_global_peak'],flush=True)
    for r in starts.values():r['gate']=gate(r)
    for r in runs.values():
        W=np.array([v['Wmax'] for v in r['rows']]);tt=np.array([v['t'] for v in r['rows']]);growth=np.gradient(np.log(W),tt)
        for i,row in enumerate(r['rows']):row['logW']=float(np.log(W[i]));row['inverseW']=float(1/W[i]);row['log_growth_rate']=float(growth[i])
        r['sampled_local_maximum_indices']=[i for i in range(1,len(W)-1) if W[i]>=W[i-1] and W[i]>W[i+1]]
    common=min(r['rows'][-1]['t'] for r in runs.values());t=np.linspace(0,common,121);comparisons={}
    for a,b in ((64,80),(80,128),(64,128)):
        ra,rb=runs[str(a)]['rows'],runs[str(b)]['rows'];av=np.interp(t,[x['t'] for x in ra],[x['Wmax'] for x in ra]);bv=np.interp(t,[x['t'] for x in rb],[x['Wmax'] for x in rb]);err=float(np.linalg.norm(av-bv)/np.linalg.norm(bv))
        comparisons[f'{a}-{b}']=dict(common_end=common,W_curve_relative_error=err,passes_5_percent=err<=.05)
    out=dict(status='initial_resolution_gate_failed',initial=[starts[n] for n in sorted(starts)],runs=runs,common_available_end=common,
             common_resolved_interval=None,comparison= comparisons,periodic_join_finding=initial['periodic_join_finding'],
             raw_curl_join_jump=initial['raw_strain_curl_trace_jump_at_x_faces_y0_z1'],new_time_evolution=False,
             decision='Do not extend time under the adaptive rule: original global peak fails the width gate already at initialization, including256.',
             data_limit='Only saved-time maxima are observable. No within-step event history exists for these archived fields.')
    save(ROOT/'results.json',out);save(ROOT/'field-manifest.json',manifest);save(ROOT/'slices.json',snapshots)
    print(out['decision'],flush=True)
if __name__=='__main__':main()
