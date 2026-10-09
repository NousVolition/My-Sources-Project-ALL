"""Postprocess saved results. No fluid or marker integration is performed."""
from pathlib import Path
import json, hashlib
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parent
def read(n):return json.loads((ROOT/n).read_text(encoding='utf-8'))
def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def wrap(x):return (x+3)%6-3
def rms(x,axis=None):return np.sqrt(np.mean(x*x,axis=axis))
def distance_matrix(x):
    d=np.linalg.norm(wrap(x[:,None]-x[None]),axis=-1);np.fill_diagonal(d,np.inf);return d
def neighbors(d):return d<=np.sort(d,axis=1)[:,7,None]+1e-12
def radius(x):return np.sqrt(np.mean(np.sum((x-x.mean(-2,keepdims=True))**2,axis=-1),axis=-1))
def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    assert read('state.json')['status']=='complete'
    meta=read('results.json');cases=meta['cases'];con=np.load(ROOT/'construction.npz',allow_pickle=False)
    checks={'source_hashes':{},'field_hashes':{},'result_hashes':{}}
    for path,h in meta['source_hashes'].items():checks['source_hashes'][path]=digest(Path(path))==h
    for path,h in meta['fields'].items():checks['field_hashes'][path]=digest(ROOT.parent/'matched-stretch-study/runs'/path)==h
    data={}
    for name,run in meta['runs'].items():
        checks['result_hashes'][name]=digest(ROOT/(name+'.npz'))==run['sha256']
        d=np.load(ROOT/(name+'.npz'),allow_pickle=False);data[name]={key:d[key] for key in d.files}
        assert np.isfinite(d['positions']).all() and np.isfinite(d['states']).all()
        assert d['positions'].shape[1:]==(90,192,3) and d['states'].shape[1:]==(90,192)
        assert d['time'][0]==0 and d['time'][-1]==.1
    assert all(all(v.values()) for v in checks.values())
    geometry={}
    for name in ['menger','koch','mixed','circle','ball']:
        x=con[name];D=distance_matrix(x);adj=neighbors(D);adj=adj|adj.T
        eig=np.linalg.eigvalsh((x-x.mean(0)).T@(x-x.mean(0))/192)
        geometry[name]={'count':len(x),'radius':float(radius(x)),'centroid':x.mean(0),'minimum_separation':D.min(),
            'components':connected_components(sparse.csr_matrix(adj),directed=False)[0],
            'degree_range':[adj.sum(1).min(),adj.sum(1).max()],'embedding_covariance_eigenvalues':eig}
        assert len(np.unique(x,axis=0))==192 and np.linalg.norm(x.mean(0))<1e-14 and abs(radius(x)-.15)<1e-14
    metrics={};symmetry={}
    for name,d in data.items():
        xs=d['positions'];zs=d['states'];rows=[];passive_error=0.;sign_error=0.
        for ci,c in enumerate(cases):
            x=xs[:,ci];z=zs[:,ci];D0=distance_matrix(x[0]);initial=neighbors(D0);DF=distance_matrix(x[-1]);final=neighbors(DF)
            retained=(initial&final).sum(1)/initial.sum(1)
            rows.append({**c,'state_rms_initial':rms(z[0]),'state_rms_final':rms(z[-1]),'state_rms_ratio':rms(z[-1])/rms(z[0]),
                'state_max_abs':abs(z).max(),'radius_initial':radius(x[0]),'radius_final':radius(x[-1]),
                'minimum_separation_final':DF.min(),'minimum_separation_sampled':min(distance_matrix(q).min() for q in x),
                'initial_neighbor_retention_final':retained.mean(),'radius_curve':radius(x),'state_rms_curve':rms(z,axis=1)})
            if c['mu']==0:
                same=next(i for i,e in enumerate(cases) if e['geometry']==c['geometry'] and e['mu']==0)
                passive_error=max(passive_error,float(abs(x-xs[:,same]).max()))
                if c['initial_state']=='positive':
                    ni=next(i for i,e in enumerate(cases) if e['geometry']==c['geometry'] and e['profile']==c['profile'] and e['initial_state']=='negative' and e['mu']==0)
                    sign_error=max(sign_error,float(abs(z+zs[:,ni]).max()))
        assert passive_error==0 and sign_error<1e-12
        symmetry[name]={'passive_positions_equal_across_relationship_conditions':passive_error,'passive_positive_negative_sign_error':sign_error}
        pairs=[]
        for ci in range(0,90,2):
            p=rows[ci];a=rows[ci+1];delta=wrap(xs[:,ci+1]-xs[:,ci]);norm=np.linalg.norm(delta,axis=-1)
            pairs.append({'geometry':p['geometry'],'profile':p['profile'],'initial_state':p['initial_state'],
                'passive_case':ci,'active_case':ci+1,'position_effect_rms_final':rms(norm[-1]),'position_effect_max_final':norm[-1].max(),
                'position_effect_rms_curve':rms(norm,axis=1),'radius_change_percent':100*(a['radius_final']/p['radius_final']-1),
                'state_difference_rms_final':rms(zs[-1,ci+1]-zs[-1,ci]),
                'neighbor_retention_difference':a['initial_neighbor_retention_final']-p['initial_neighbor_retention_final']})
        metrics[name]={'times':d['time'],'cases':rows,'pairs':pairs}
    sensitivity={}
    for label,refname,othername in [('step','step64','base64'),('grid','step64','grid128'),('snapshot','step64','snapshot_coarse')]:
        ref=data[refname];other=data[othername];idx=[int(np.argmin(abs(ref['time']-t))) for t in other['time']]
        assert np.allclose(ref['time'][idx],other['time'],rtol=0,atol=1e-12)
        norm=np.linalg.norm(wrap(other['positions']-ref['positions'][idx]),axis=-1)
        dz=rms(other['states']-ref['states'][idx],axis=2)/np.maximum(rms(ref['states'][idx],axis=2),1e-14)
        rows=[{'id':c['id'],'max_position_difference':norm[:,i].max(),'final_position_rms_difference':rms(norm[-1,i]),
               'max_relative_state_rms_difference':dz[:,i].max()} for i,c in enumerate(cases)]
        sensitivity[label]={'reference':refname,'comparison':othername,'cases':rows,'max_position_difference':norm.max(),
            'max_relative_state_rms_difference':dz.max(),'max_final_position_rms_difference':rms(norm[-1],axis=1).max()}
        if label=='step':
            for r in rows:r['pass']=bool(r['max_position_difference']<=.001 and r['max_relative_state_rms_difference']<=.005)
            sensitivity[label]['passes']=sum(r['pass'] for r in rows)
        rp=metrics[refname]['pairs'];op=metrics[othername]['pairs']
        sensitivity[label]['radius_effect_sign_agreement']=sum(np.sign(a['radius_change_percent'])==np.sign(b['radius_change_percent']) for a,b in zip(rp,op))
        sensitivity[label]['max_radius_effect_difference_percentage_points']=max(abs(a['radius_change_percent']-b['radius_change_percent']) for a,b in zip(rp,op))
    checks['symmetry']=symmetry
    # Independent vectorized matrix assembly for one full network RHS.
    from study import relationship_rhs
    ci=47;x=data['step64']['positions'][-1,ci];z=data['step64']['states'][-1,ci]
    all_edges=con['edges'];ed=all_edges[(all_edges[:,0]//192)==ci]-192*ci;i,j=ed.T
    rel=wrap(x[j]-x[i]);w=np.exp(-np.minimum(np.sum(rel*rel,axis=1)/.02,50));W=np.zeros((192,192));W[i,j]=w;W[j,i]=w
    A=con['a'][ci*192:(ci+1)*192];B=con['b'][ci*192:(ci+1)*192];M=con['mobility'][ci*192:(ci+1)*192]
    dz=A*z+B*(W@z)/W.sum(1)
    diff=wrap(x[None]-x[:,None]);factor=W*np.tanh((z[:,None]+z[None])/2)
    drift=M[:,None]*np.sum(factor[:,:,None]*diff,axis=1)
    got,gz=relationship_rhs(x,z,A,B,M,ed)
    checks['independent_rhs_max_error']=max(abs(got-drift).max(),abs(gz-dz).max());assert checks['independent_rhs_max_error']<1e-12
    checks['finite_arrays']=True;checks['all_expected_cases_and_times']=True
    output={'geometry':geometry,'metrics':metrics,'sensitivity':sensitivity,'checks':checks,'analysis_source_sha256':digest(Path(__file__))}
    (ROOT/'analysis.json').write_text(json.dumps(output,indent=2,default=plain,allow_nan=False)+'\n',encoding='utf-8')
    brief={'status':'complete','cases':360,'matched_pairs_per_configuration':45,'step_screen_passes':sensitivity['step']['passes'],
        'step_max_position':sensitivity['step']['max_position_difference'],'grid_max_position':sensitivity['grid']['max_position_difference'],
        'snapshot_max_position':sensitivity['snapshot']['max_position_difference'],
        'radius_effect_sign_agreement_grid':sensitivity['grid']['radius_effect_sign_agreement'],
        'effects_step64':metrics['step64']['pairs']}
    (ROOT/'summary.json').write_text(json.dumps(brief,indent=2,default=plain,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in brief.items() if k!='effects_step64'},default=plain))

if __name__=='__main__':main()
