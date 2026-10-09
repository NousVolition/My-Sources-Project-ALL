"""Independent saved-array and artifact checks for this local pilot."""
from pathlib import Path
from html.parser import HTMLParser
import json,hashlib
import numpy as np
from scipy import fft,ndimage
ROOT=Path(__file__).resolve().parent
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(n):return json.loads((ROOT/n).read_text(encoding='utf-8'))

def main():
    result=read('results.json');analysis=read('analysis.json');cases=result['cases'];construction=np.load(ROOT/'construction.npz')
    assert len(cases)==len({c['id'] for c in cases})==90
    largest_error=0.;arrays=0
    for name in result['runs']:
        d=np.load(ROOT/(name+'.npz'));x=d['positions'];z=d['states'];assert np.array_equal(z[0],construction['initial_states'])
        assert np.isfinite(x).all() and np.isfinite(z).all();arrays+=3
        for j,c in enumerate(cases):
            assert np.array_equal(x[0,j],construction[c['geometry']])
        for p in analysis['metrics'][name]['pairs']:
            ip,ia=p['passive_case'],p['active_case'];cp,ca=cases[ip],cases[ia]
            assert cp['mu']==0 and ca['mu']==20
            assert all(cp[k]==ca[k] for k in ['geometry','profile','initial_state'])
            def size(points):
                # Pair-distance identity provides an independent radius formula.
                delta=points[:,None]-points[None,:]
                return np.sqrt(np.sum(delta*delta)/(2*192**2))
            effect=100*(size(x[-1,ia])/size(x[-1,ip])-1)
            largest_error=max(largest_error,abs(effect-p['radius_change_percent']))
    assert largest_error<1e-10
    from study import interp,relationship_rhs,graph,geometries
    pos=geometries()['mixed'];edges=graph(pos);u=np.array([1.,-2.,.5]);field=np.broadcast_to(u[:,None,None,None],(3,8,8,8))
    drift,dz=relationship_rhs(pos,np.zeros(192),np.full(192,-20.),np.full(192,10.),np.ones(192),edges)
    vel=interp(field,pos)+drift;shifted=pos+.1*vel
    translation_error=float(abs(shifted-pos-.1*u).max());assert translation_error<1e-12 and not np.any(dz)
    h=np.load(ROOT.parent/'matched-stretch-study/runs/baseline-n64-base/field-000.npy');field=fft.irfftn(h,s=(64,64,64),axes=(-3,-2,-1),workers=1)
    independent=np.stack([ndimage.map_coordinates(c,(((pos+3)%6)*64/6).T,order=1,mode='grid-wrap') for c in field],axis=1)
    interp_error=float(abs(interp(field,pos)-independent).max());assert interp_error<1e-11
    # Any row-normalized graph preserves a uniform neighbor state.
    _,got=relationship_rhs(pos,np.full(192,.2),np.full(192,-20.),np.full(192,30.),np.zeros(192),edges)
    uniform_error=float(abs(got-2).max());assert uniform_error<1e-12
    class Links(HTMLParser):
        def __init__(self):super().__init__();self.paths=[]
        def handle_starttag(self,tag,attrs):
            for key,value in attrs:
                if key in ['src','href'] and not value.startswith(('http','#')):self.paths.append(value)
    parser=Links();parser.feed((ROOT/'report.html').read_text(encoding='utf-8'))
    for path in parser.paths:
        if path!='verification.json':assert (ROOT/path).is_file(),path
    sign_flips=[]
    for p,q in zip(analysis['metrics']['step64']['pairs'],analysis['metrics']['grid128']['pairs']):
        if np.sign(p['radius_change_percent'])!=np.sign(q['radius_change_percent']):
            sign_flips.append({'geometry':p['geometry'],'profile':p['profile'],'initial_state':p['initial_state'],
                'radius_effect_64_percent':p['radius_change_percent'],'radius_effect_128_percent':q['radius_change_percent']})
    artifact_hashes={p.name:sha(p) for p in ROOT.iterdir() if p.is_file() and p.name not in ['verification.json','verification.tmp']}
    record={'status':'passed','trajectory_arrays_checked':arrays,'cases_per_configuration':90,'configurations':4,
        'independent_radius_effect_max_error_percentage_points':largest_error,'uniform_translation_error':translation_error,
        'independent_saved_velocity_interpolation_error':interp_error,'uniform_network_state_rhs_error':uniform_error,
        'all_report_local_links_exist':True,'source_and_field_hashes_verified_by_analysis':True,
        'source_hash_check_count':len(analysis['checks']['source_hashes']),'field_hash_check_count':len(analysis['checks']['field_hashes']),
        'grid_radius_effect_sign_flips':sign_flips,'artifact_sha256':artifact_hashes,
        'corrected_postprocessing_issue':'The first analysis invocation saved correct JSON outputs but its final console summary failed on NumPy int64 serialization. Console serialization was fixed; analysis then completed with exit code 0. No trajectory reruns or numerical changes.',
        'numerical_limitations':'Tracking step screen passed. Grid and snapshot sensitivities remain. Original nonsmooth periodic initial fluid is unchanged; no physical accuracy claim.'}
    (ROOT/'verification.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ['artifact_sha256','grid_radius_effect_sign_flips']}))

if __name__=='__main__':main()
