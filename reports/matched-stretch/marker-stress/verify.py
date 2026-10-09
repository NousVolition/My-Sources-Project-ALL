"""Independent saved-trajectory and endpoint-velocity checks; no integration."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.ndimage import map_coordinates
ROOT=Path(__file__).resolve().parent;STUDY=ROOT.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def wrap(x):return (x+3)%6-3
def main():
    data=json.loads((ROOT/'results.json').read_text());assert len(data['fields'])==93 and len(data['runs'])==8
    assert all(sha(STUDY/n)==digest for n,digest in data['source_hashes'].items())
    assert json.loads((STUDY/'protocol.json').read_text())['external_force']==0
    assert sha(ROOT/'run.py')==data['run_script_sha256'] and sha(ROOT/'refine.py')==data['refine_script_sha256']
    endfields={};checks=[];ij=np.array(data['edges'])
    for name,record in data['runs'].items():
        path=ROOT/(name+'.npz');assert sha(path)==record['trajectory_sha256'];z=np.load(path,allow_pickle=False);assert np.isfinite(z['positions']).all() and np.isfinite(z['velocities']).all()
        assert np.array_equal(z['positions'][:,500:],z['positions'][:,:5]);rid=record['source'];n=record['n']
        if rid not in endfields:
            p=STUDY/'runs'/rid/'field-030.npy';assert sha(p)==data['fields'][rid+'/field-030.npy'];h=np.load(p,allow_pickle=False);endfields[rid]=np.fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1))
        pos=z['positions'][-1];coordinates=(((pos+3)%6)*n/6).T;sample=np.stack([map_coordinates(u,coordinates,order=1,mode='grid-wrap') for u in endfields[rid]],axis=1)
        error=float(np.max(abs(sample-z['velocities'][-1])));assert error<1e-10
        max_metric_error=0.
        for k,t in enumerate(z['time']):
            for g in (0,1):
                p=z['positions'][k];v=z['velocities'][k];subset=slice(g*125,(g+1)*125);edge=ij[g*300:(g+1)*300];mean=np.sum(v[subset],axis=0)/125;variance=sum(float(np.dot(q-mean,q-mean)) for q in v[subset])/125
                d=np.sqrt(np.sum(wrap(p[edge[:,1]]-p[edge[:,0]])**2,axis=1));delta=wrap(p[250+g*125:250+(g+1)*125]-p[subset]);amp=np.sqrt(np.mean(np.sum(delta*delta,axis=1)))/.0009375
                r=record['metrics'][k][g];e=max(abs(variance**.5-r['velocity_deviation_rms']),abs(float(np.median(d))-r['neighbor_distance_quantiles'][2]),abs(amp-r['seed_amplification_rms']));assert e<1e-9;max_metric_error=max(max_metric_error,e)
        checks.append({'name':name,'snapshots':len(z['time']),'trajectory_sha256':sha(path),'endpoint_independent_interpolation_max_error':error,'recomputed_group_metrics_max_error':max_metric_error,'identical_clones_agree':True})
    verification={'status':'passed','scope':'Saved-track, diagnostic and endpoint sampling consistency, not physical convergence','saved_source_fields_hashed':93,'tracking_configurations':8,'main_labels':250,'position_perturbations':250,'identical_clones':5,'new_fluid_simulations':0,'unchanged_numerical_fluid_sources':data['source_hashes'],'analytic_checks':data['analytic_checks'],'runs':checks}
    (ROOT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n');print(json.dumps({'status':'passed','maximum_independent_velocity_error':max(x['endpoint_independent_interpolation_max_error'] for x in checks),'maximum_metric_error':max(x['recomputed_group_metrics_max_error'] for x in checks)}))
if __name__=='__main__':main()
