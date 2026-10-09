"""Audit existing starts and first steps. No long simulation or source mutation.

Use --study-root to locate the saved matched-stretch-study directory.
"""
from pathlib import Path
import argparse, ast, hashlib, json, os, sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source'))
import clay_hug as m
from diagnostics import Meter, initial, mirror, norm, dot
from numerics import Flow
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,x):(ROOT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def measure(u,ops,phi):
    u=np.asarray(u);h=np.fft.fftn(u,axes=(-3,-2,-1));k=np.asarray(ops[:3]);wh=1j*np.stack((k[1]*h[2]-k[2]*h[1],k[2]*h[0]-k[0]*h[2],k[0]*h[1]-k[1]*h[0]));w=np.fft.ifftn(wh,axes=(-3,-2,-1)).real
    # Raw component diagnostics deliberately use every represented mode.
    div=np.fft.ifftn(1j*np.sum(k*h,axis=0)).real;a=u-mirror(u);U=norm(u)
    return dict(E=float(norm(a)/U) if U else 0.,D=dot(a,phi) if phi is not None else None,U=float(U),A=float(norm(a)),W=float(np.sqrt(np.sum(w*w,axis=0)).max()),energy=.5*6**3*dot(u,u),enstrophy=.5*6**3*dot(w,w),divergence_max=float(abs(div).max()),mean=u.mean(axis=(1,2,3)).tolist())
def slopes(rows,end=.04):
    a=[x for x in rows if x['t']<=end+1e-12]
    if len(a)<2:a=rows[:2]
    t=np.array([x['t'] for x in a]);out=dict(window=[float(t[0]),float(t[-1])],samples=len(a),kind='Least-squares slope over saved samples; not an instantaneous derivative')
    for k in ('E','D','W'):
        if all(x.get(k) is not None for x in a):out[k]=float(np.polyfit(t,[x[k] for x in a],1)[0])
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--study-root',type=Path,default=ROOT.parent);args=ap.parse_args();study=args.study_root.resolve()
    assert sha(ROOT/'source/clay_hug.py')=='422b7160f9e96badfc68fb76958c3cb512982905b73b37f8db92e5fd37379c01'
    assert sha(ROOT/'source/numerics.py')=='22b169971aaf2efd56d47f685495a575e1622edd56c199b07a785fc97e952f46'
    components={};starts={};steps={};histories={};checks={};inputs={};pair_checks={}
    (ROOT/'raw-fields').mkdir(exist_ok=True)
    def track(path):inputs[str(path.relative_to(study))]=sha(path)
    def component(name,u,ops,phi):
        p=ROOT/'raw-fields'/f'{name}.npz';np.savez_compressed(p,u=np.asarray(u));components[name]=dict(measurements=measure(u,ops,phi),file=str(p.relative_to(ROOT)),sha256=sha(p))
    def start_and_step(name,u,ops,phi,dt,pre=None):
        u=np.asarray(u);v=np.asarray(m.project(u,*ops[:5]));s=np.asarray(m.advance(v,dt,ops));starts[name]=dict(before_projection=measure(u if pre is None else pre,ops,phi),after_projection=measure(v,ops,phi),projection_relative_change=float(norm(v-u)/norm(u)))
        first=measure(s,ops,phi);first.update(t=dt);steps[name]=dict(dt=dt,force_kind='none',force_norm_at_time_zero=0.,first=first,first_difference_from_unforced=0.)
        return v,s
    # Exact mirrored source construction, with the unprojected components exposed.
    b,phi,ops,c=initial(m,33);x,y,z,*_=m.grids(33);psi=np.exp(-12*((x-.9)**2+y*y))*(x>0)
    raw_a=np.asarray([m.deriv(psi,ops[1],ops[4]),-m.deriv(psi,ops[0],ops[4]),np.zeros_like(psi)])
    raw_ab=(raw_a+mirror(raw_a))*.5;b=(b+mirror(b))*.5
    raw_c=np.zeros_like(b);raw_c[1]=.25*np.exp(-8*((x-.55)**2+(y-.4)**2+z*z));d=mirror(c)
    for name,u in [('one_half_filtered_curl',raw_a),('AB_before_projection',raw_ab),('C_before_projection',raw_c),('C_after_projection',c),('D_after_projection',d),('odd_template',phi)]:component(name,u,ops,phi)
    reference_path=study/'received-mirror-abc/complete-fluid-scripts/runs/fluid_sources/fluid-sources.json';reference=read(reference_path);track(reference_path)
    sourcefields={'ab':b,'ab_cd':b+.5*c+.5*d,'ab_c':b+c,'ab_uneven':b+.7*c+.3*d,'ab_c_mirror':b+d}
    source_one={}
    for label,u in sourcefields.items():
        dt=.5/np.ceil(.5/(.05*ops[-1]/np.sqrt(np.sum(u*u,axis=0)).max()))
        v,s=start_and_step(label,u,ops,phi,float(dt));source_one[label]=s;component(label+'_before_final_projection',u,ops,phi)
        r=reference[label];errs={k:abs(starts[label]['after_projection'][k]-r[0][k]) for k in ('E','D','W')};assert max(errs.values())<1e-10
        checks[label]=dict(initial_reference_errors=errs);histories[label]=dict(rows=r,slopes=slopes(r),sampling_note='Original times rounded to 0.001; first stored interval is about 0.1, not 0.04.')
    pair_checks['C_vs_reflected_C']=dict(first_field_error=float(norm(source_one['ab_c_mirror']-mirror(source_one['ab_c']))/norm(source_one['ab_c'])),max_E_difference=max(abs(a['E']-b['E']) for a,b in zip(reference['ab_c'],reference['ab_c_mirror'])),max_D_sum=max(abs(a['D']+b['D']) for a,b in zip(reference['ab_c'],reference['ab_c_mirror'])))
    # Existing unforced controls: zero bias and opposite signed starts, all three grids.
    for n in (33,49,65):
        b,p,ops,c=initial(m,n);one={}
        for sign in (0,1,-1):
            rid=f'n{n}-s{sign}-shape0-none-base';rp=study/'methods-controls/runs'/rid/'result.json';r=read(rp);track(rp);assert r['status']=='complete'
            u=b+sign*.01*norm(b)*p;label=f'mirror-{n}-{sign}';v,s=start_and_step(label,u,ops,p,r['dt']);one[sign]=s
            f=rp.parent/r['fields'][0]['file'];track(f);saved=np.load(f)['u'];err=norm(saved-u)/norm(u);assert err<1e-12
            checks[label]=dict(initial_saved_field_relative_error=float(err))
            rows=[dict(t=v['t'],**{k:v[k] for k in ('E','D','W')}) for v in r['rows']];histories[label]=dict(rows=rows,slopes=slopes(rows))
        pair_checks[f'mirror-{n}']=dict(first_field_error=float(norm(one[-1]-mirror(one[1]))/norm(one[1])),first_E_difference=measure(one[1],ops,p)['E']-measure(one[-1],ops,p)['E'],first_D_sum=measure(one[1],ops,p)['D']+measure(one[-1],ops,p)['D'])
    # Reconstruct the exact older breathing setup without executing its run loop.
    bp=study/'breathing-rerun-intake/v2/reproduction/breathing_rerun.py';track(bp);tree=ast.parse(bp.read_text());stop=next(i for i,v in enumerate(tree.body) if isinstance(v,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='out' for t in v.targets));tree.body=tree.body[:stop]
    ns={'__file__':str(bp),'__name__':'breathing_initialization_only'};exec(compile(tree,str(bp),'exec'),ns)
    ops=ns['ops'];p=np.asarray(ns['phi']);force=np.asarray(ns['breath_u']);base=np.asarray(ns['base'])
    raw_force=np.stack((np.exp(-8*(ns['x']**2+ns['y']**2))*np.sign(ns['x']),np.zeros_like(ns['x']),np.zeros_like(ns['x'])))
    component('old_breathing_before_projection',raw_force,ops,p);component('old_breathing_after_projection',force,ops,p)
    for label,u in ns['cases'].items():
        folder=study/'breathing-rerun-intake/v2/timeline'/label.replace(' ','-');rp=folder/'result.json';r=read(rp);track(rp);manifest=read(folder/'field-manifest.json');dt=r['dt'];key='old-breathing-'+label
        u=np.asarray(m.project(np.asarray(u),*ops[:5]));unforced=np.asarray(m.advance(u,dt,ops));coefficient=.15*np.sin(2*np.pi*dt/.4);kicked=u+coefficient*force;v=np.asarray(m.advance(kicked,dt,ops))
        rows=[]
        for i,item in enumerate(manifest):
            if r['series'][i]['t']>.04 and i>1:break
            fp=folder/item['file'];assert sha(fp)==item['sha256'];track(fp);a=np.load(fp)['u'];obs=measure(a,ops,p);obs['t']=r['series'][i]['t'];rows.append(obs)
            if i==1:err=float(norm(a-v)/norm(v));assert err<1e-12;checks[key]=dict(first_saved_field_relative_error=err)
        starts[key]=dict(after_projection=measure(u,ops,p));steps[key]=dict(dt=dt,force_kind='Per-step velocity kick, not multiplied by dt',force_norm_at_time_zero=0.,first_kick_time=dt,first_kick_coefficient=float(coefficient),first_kick_norm=float(norm(coefficient*force)),kick_projection_onto_D=float(dot(coefficient*(force-mirror(force)),p)),first_difference_from_unforced=float(norm(v-unforced)/norm(unforced)),first=dict(t=dt,**measure(v,ops,p)))
        histories[key]=dict(rows=rows,slopes=slopes(rows));
    # New driven experiments use both Heun stages. Show that F(0)=0 does not mean no force in the first step.
    for parity in ('even-gap','odd-lean'):
        rid=f'n49-{parity}-A0.5-P0.2-s1-base';rp=study/'model-layers/fluid-runs'/rid/'result.json';r=read(rp);track(rp);assert r['status']=='complete'
        b,p,ops,c=initial(m,49);u=b+.01*norm(b)*p;x,y,z,*_=m.grids(49);raw=np.zeros_like(u);raw[0]=np.exp(-8*(x*x+y*y))*np.sign(x)
        gap=np.asarray(m.project(raw,*ops[:5]));gap=(gap+mirror(gap))*.5;gap-=gap.mean(axis=(1,2,3),keepdims=True);gap/=norm(gap)
        shape=gap if parity=='even-gap' else p.copy();shape-=shape.mean(axis=(1,2,3),keepdims=True);shape=np.asarray(m.project(shape,*ops[:5]));shape/=norm(shape);F=.5*norm(b)/.4*shape;dt=r['dt'];F1=np.zeros_like(F);F2=np.sin(2*np.pi*dt/.2)*F
        rhs1=np.asarray(m.rhs(u,*ops[:5]))+F1;y=u+dt*rhs1;rhs2=np.asarray(m.rhs(y,*ops[:5]))+F2;v=np.asarray(m.project(u+.5*dt*(rhs1+rhs2),*ops[:5]));unforced=np.asarray(m.advance(u,dt,ops));key='new-drive-'+parity
        obs=measure(v,ops,p);saved=r['dense'][1];err=max(abs(obs[k]-saved[k])/max(1,abs(saved[k])) for k in ('E','D','W','energy','enstrophy'));assert err<1e-11
        checks[key]=dict(first_dense_record_scaled_error=err);starts[key]=dict(after_projection=measure(u,ops,p))
        steps[key]=dict(dt=dt,force_kind='Acceleration integrated at both Heun stages',force_norm_at_time_zero=0.,force_norm_second_stage=float(norm(F2)),velocity_increment_from_force_norm=float(norm(.5*dt*F2)),direct_D_force_at_second_stage=dot(F2-mirror(F2),p),first_difference_from_unforced=float(norm(v-unforced)/norm(unforced)),first=dict(t=dt,**obs))
        histories[key]=dict(rows=r['dense'],slopes=slopes(r['dense']))
    # Matched stretch is a different field. Reconstruct raw components and its first step at 64.
    f=Flow(64,.001,workers=2);axis=np.arange(64)*f.dx-3;x,y,z=np.meshgrid(axis,axis,axis,indexing='ij',sparse=True);env=np.exp(-.04*(x*x+y*y+z*z));raw=np.asarray([np.broadcast_to(-40*x*env,(64,)*3),np.broadcast_to(-40*y*env,(64,)*3),np.broadcast_to(80*z*env,(64,)*3)])
    psi=np.broadcast_to(.4*np.exp(-(x*x+y*y)/.2**2),(64,)*3);ph=f.hat(psi)*f.keep;tube=np.asarray([f.real(1j*f.k[1]*ph),f.real(-1j*f.k[0]*ph),np.zeros((64,)*3)])
    gx,gy,gz,kx,ky,kz,k2,keep,dx=m.grids(64);ops=(kx,ky,kz,k2,keep,dx)
    component('matched_raw_strain',raw,ops,None);component('matched_filtered_tube',tube,ops,None);component('matched_before_projection',raw+tube,ops,None)
    rp=study/'runs/baseline-n64-base/result.json';r=read(rp);track(rp);h=f.project(f.hat(raw+tube));saved=np.load(rp.parent/'field-000.npy');err=np.linalg.norm(h-saved)/np.linalg.norm(saved);assert err<1e-12;checks['matched-64']=dict(initial_spectral_field_error=float(err))
    v,accum,speed=f.step(h,r['dt']);obs=f.observe(v,r['dt']);starts['matched-64']=dict(before_projection=components['matched_before_projection']['measurements'],after_projection=measure(f.real(h),ops,None));steps['matched-64']=dict(dt=r['dt'],force_kind='none',force_norm_at_time_zero=0.,first=dict(t=r['dt'],W=obs['Wmax'],energy=obs['energy'],enstrophy=obs['enstrophy'],divergence_max=obs['divergence_max']))
    histories['matched-64']=dict(rows=[dict(t=v['t'],W=v['Wmax']) for v in r['series']],slopes=slopes([dict(t=v['t'],W=v['Wmax']) for v in r['series']]))
    resolution={}
    for key,p in [('mirror',study/'methods-controls/analysis.json'),('matched',study/'adaptive-peak-audit/results.json')]:
        a=read(p);track(p)
        resolution[key]={k:a[k] for k in (('curve_comparisons','mirror_pairs') if key=='mirror' else ('initial','comparison','common_resolved_interval','decision'))}
    assert max(x['first_field_error'] for x in pair_checks.values())<1e-12
    save('results.json',dict(status='complete',scope='Existing starts, reconstructed first steps, saved early histories and completed refinement audits. No new long evolution.',definitions={'E':'norm(u-Mu)/norm(u), unsigned','D':'mean spatial inner product of u-Mu with unit antisymmetric phi, signed','divergence':'Spatial incompressibility residual; different from separation between experiments','raw_measurement':'Full represented FFT spectrum before final projection. Derivative-built input components already include their original filtering.'},components=components,starts=starts,steps=steps,histories=histories,pairs=pair_checks,resolution=resolution,verification=checks,input_hashes=inputs,source_hashes={p.name:sha(p) for p in (ROOT/'source').glob('*.py')}))
    print(json.dumps(dict(starts=len(starts),first_steps=len(steps),components=len(components),paired_checks=len(pair_checks),source_files=len(inputs)),indent=2),flush=True)
if __name__=='__main__':main()
