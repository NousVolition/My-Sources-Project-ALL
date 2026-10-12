"""Intervention screening from the completed 64-cubed stationary candidate."""
import argparse,time
import numpy as np
from solver import Flow,resize
from simulate import response
from common import ROOT,protocol,digest,sources,save,now

def main(label):
 p=protocol();source=ROOT/'data/screen/n64_nu0.025_dt0.025_seed300000.npz';c=p.copy();c['n']=p['fine_n'] if label=='fine' else p['n'];c['response_dt']=p['response_dt'] if label=='base' else p['response_dt']/2
 with np.load(source) as z:old=z['final']
 f=Flow(c['n'],c['nu'],c['power']);h=f.project(resize(old,c['n']));start=time.time();started=now();z=response(f,h,c);out=ROOT/'data/screen'/f'response_{label}_300000.npz';np.savez_compressed(out,**z)
 q=z['response_checks'];g=p['gates'];checks=dict(cfl=float(z['response_cfl'])<=g['cfl_max'],resolution=float(q[:,0].min())>=g['keta_min'],tail=float(q[:,1].max())<=g['tail_energy_max'],total_budget=float(q[:,2].max())<=g['budget_relative_max'],difference_budget=float(q[:,3].max())<=g['budget_relative_max'])
 save(str(out.with_suffix('.json').relative_to(ROOT)),dict(label=label,config={k:c[k] for k in ['n','nu','power','response_dt','horizon','strengths','impulse_width']},input_path=source.relative_to(ROOT).as_posix(),input_sha256=digest(source),source_sha256={**sources(),'screen_response.py':digest(ROOT/'screen_response.py')},started_utc=started,completed_utc=now(),elapsed_seconds=time.time()-start,sha256=digest(out),checks=checks,minimum_keta=float(q[:,0].min()),maximum_cfl=float(z['response_cfl'])))
 print(label,checks,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('label',choices=['base','half','fine']);a=p.parse_args();main(a.label)
