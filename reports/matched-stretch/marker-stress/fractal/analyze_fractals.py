"""Count boxes at finite scales and compare passive-label reconstructions."""
from pathlib import Path
import json,hashlib,os
import numpy as np
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def wrap(x):return (x+3)%6-3
def count(p,epsilon,offset=0.):return len(np.unique(np.floor(p/epsilon+offset).astype(np.int64),axis=0))
def profile(p,base):
    # Remove translation, keep the fixed original scale; do not refit size.
    q=(p-p.mean(axis=0))/.375;counts=[]
    for offset in (0.,.25,.5,.75):counts.append([count(q,base**(-float(k)),offset) for k in range(7)])
    c=np.array(counts);local=np.diff(np.log(c),axis=1)/np.log(base);fits={}
    for a,b in ((1,3),(2,4)):
        k=np.arange(a,b+1);fit=[]
        for row in c:
            slope,intercept=np.polyfit(k*np.log(base),np.log(row[a:b+1]),1);pred=slope*k*np.log(base)+intercept;ss=np.sum((np.log(row[a:b+1])-np.mean(np.log(row[a:b+1])))**2);r2=1-np.sum((pred-np.log(row[a:b+1]))**2)/ss if ss>1e-25 else None
            fit.append({'slope':float(slope),'r_squared':None if r2 is None else float(r2)})
        fits[f'{a}-{b}']=fit
    return {'base':base,'k':list(range(7)),'count':c.tolist(),'local_slopes':local.tolist(),'fits':fits,'marker_count':len(p),'finest_count_fraction':(c[:,-1]/len(p)).tolist()}
def main():
    datasets=[];f=json.loads((ROOT/'results.json').read_text());o=json.loads((ROOT/'menger-results.json').read_text());qdata=json.loads((ROOT/'quarter-results.json').read_text())
    datasets.append(('carpet',f,[{'name':'carpet3','start':0,'stop':512,'base':3,'depth':3},{'name':'carpet4','start':512,'stop':4608,'base':3,'depth':4}],''))
    datasets.append(('menger',o,o['groups'],'menger-'));datasets.append(('quarter',qdata,qdata['groups'],'quarter-'));out={'profiles':{},'comparisons':{},'references':{},'verification':{'status':'passed','new_fluid_simulations':0,'raw_files':{}}}
    for label,data,groups,prefix in datasets:
        raw={}
        for name,r in data['runs'].items():
            p=ROOT/(prefix+name+'.npz');assert sha(p)==r['trajectory_sha256'];z=np.load(p,allow_pickle=False);assert np.isfinite(z['positions']).all() and z['time'][-1]==.1;raw[name]=z;out['verification']['raw_files'][p.name]=sha(p)
        for group in groups:
            key=group['name'];sl=slice(group['start'],group['stop']);out['profiles'][key]={}
            for name,z in raw.items():
                assert np.array_equal(z['positions'][0],raw['n64']['positions'][0]);out['profiles'][key][name]={}
                for t in (0.,.1):
                    p=z['positions'][np.argmin(abs(z['time']-t)),sl];out['profiles'][key][name][str(t)]={'3d':profile(p,group['base']),'xy':profile(p[:,:2],group['base'])}
            out['comparisons'][key]={}
            for title,a,b in [('tracking_step','n64','n64_half'),('fluid_grid','n64_half','n128'),('snapshot_spacing','n64_half','coarse_saved_times')]:
                za,zb=raw[a],raw[b];ia=np.array([np.argmin(abs(za['time']-t)) for t in zb['time']]);assert np.max(abs(za['time'][ia]-zb['time']))<1e-12
                aa=za['positions'][ia,sl];bb=zb['positions'][:,sl];d=np.linalg.norm(wrap(aa-bb),axis=2);shape=np.linalg.norm((aa-aa.mean(axis=1,keepdims=True))-(bb-bb.mean(axis=1,keepdims=True)),axis=2)
                out['comparisons'][key][title]={'maximum_position_difference':float(d.max()),'final_rms_position_difference':float(np.sqrt(np.mean(d[-1]**2))),'final_centroid_aligned_shape_rms':float(np.sqrt(np.mean(shape[-1]**2)))}
    # Exact known construction counts are separate from the shifted finite-scale fits.
    construction=np.load(ROOT/'construction.npz');other=np.load(ROOT/'menger-construction.npz')
    cases=[('Random eight-of-nine carpet',(construction['level4']+.5)/81,3,8,4),('Menger sponge',(other['menger3']+.5)/27,3,20,3)]
    cantor=np.array([0])
    for k in range(8):cantor=(cantor[:,None]*3+np.array([0,2])).ravel()
    cases.append(('Cantor set',((cantor+.5)/3**8)[:,None],3,2,8))
    line=((np.arange(81)+.5)/81)[:,None];square=np.array(np.meshgrid(line[:,0],line[:,0],indexing='ij')).reshape(2,-1).T
    cases.extend([('Filled line',line,3,3,4),('Filled square',square,3,9,4)])
    for name,p,base,keep,depth in cases:
        counts=[count(p,base**(-float(k))) for k in range(depth+3)];assert counts[:depth+1]==[keep**k for k in range(depth+1)]
        assert counts[-1]==len(p);slope=float(np.polyfit(np.arange(depth+1)*np.log(base),np.log(counts[:depth+1]),1)[0]);expected=float(np.log(keep)/np.log(base));assert abs(slope-expected)<1e-12
        out['references'][name]={'base':base,'retained_per_parent':keep,'depth':depth,'counts':counts,'expected_limiting_dimension':expected,'construction_scale_slope':slope,'finite_sample_saturates':True}
    quarter=np.load(ROOT/'quarter-construction.npz');dimension=qdata['limiting_dimension'];pressure=[]
    for d in range(1,11):
        intervals=quarter['level'+str(d)];length=intervals[:,1]-intervals[:,0];assert len(length)==2**d and abs(length.sum()-.75**d)<1e-12
        value=float(np.sum(length**dimension));assert abs(value-1)<1e-12;pressure.append(value)
    out['quarter_reference']={'dimension':dimension,'equation':'(1/4)^d+(1/2)^d=1','interval_length_power_sums':pressure,'deleted_quarter_by_stage':qdata['deleted_quarter_by_stage'],'note':'Two unequal connected pieces per interval; dimension is log(phi)/log2, not log3/log4. Coin flips exchange left/right placement of the same two contraction ratios.'}
    out['existing_patterns']={}
    for name in ['reference','grid128_refined']:
        z=np.load(ROOT.parent/(name+'.npz'),allow_pickle=False);out['existing_patterns'][name]={}
        for g,label in enumerate(['central','periodic_join']):
            out['existing_patterns'][name][label]={}
            for t in (0.,.1,.3):
                p=z['positions'][np.argmin(abs(z['time']-t)),g*125:(g+1)*125];out['existing_patterns'][name][label][str(t)]={'3d':profile(p,3),'xy':profile(p[:,:2],3)}
    out['scope']='Finite-scale counts in centered marker clouds at fixed original length scale0.375. Four diagonal grid-origin shifts; no proof of a limiting dimension from finite trajectories.'
    out['sources']=['https://www.math.stonybrook.edu/~scott/Book331/Fractal_Dimension.html','https://link.springer.com/article/10.1007/s00209-019-02426-2','https://www.math.uwaterloo.ca/~ervrscay/courses/amath343docs/week10.pdf']
    (ROOT/'analysis.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'reference_dimensions':{n:r['construction_scale_slope'] for n,r in out['references'].items()},'comparisons':out['comparisons']},indent=2))
if __name__=='__main__':main()
