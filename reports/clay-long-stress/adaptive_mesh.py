"""Replace an under-resolved coarse-start result by a verified finer solution.

This is convergence-driven global h-refinement, not a correction factor for a
fixed coarse mesh. Every level restarts from the same initial condition; rest
shape and Lorenz history are never interpolated between different meshes.
"""
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
ROOT=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def compare(coarse,fine):
    coarse,fine=Path(coarse),Path(fine)
    a=json.loads(coarse.read_text());b=json.loads(fine.read_text())
    pa=dict(a['parameters']);pb=dict(b['parameters'])
    result=dict(coarse=coarse.name,fine=fine.name,coarse_mesh=[pa['angles'],pa['layers']],fine_mesh=[pb['angles'],pb['layers']],tolerance=.05,gap_floor=.01)
    result['refines_both_axes']=pb['angles']>pa['angles'] and pb['layers']>pa['layers']
    pa={k:v for k,v in a['parameters'].items() if k not in ['angles','layers']}
    pb={k:v for k,v in b['parameters'].items() if k not in ['angles','layers']}
    if pa!=pb or a['initial']!=b['initial']:raise ValueError('Different physical setup or initial state')
    if a.get('start_time',0)!=b.get('start_time',0) or a['requested_end']!=b['requested_end']:
        raise ValueError('Compare equal requested intervals; do not truncate a failed run into a pass')
    aa=np.load(coarse.with_suffix('.npz'));bb=np.load(fine.with_suffix('.npz'))
    if not all(np.isfinite(z['state']).all() and np.all(np.diff(z['time'])>0) for z in [aa,bb]):raise ValueError('Invalid raw recording')
    pairs=[(i,int(np.argmin(abs(bb['time']-t)))) for i,t in enumerate(aa['time']) if np.min(abs(bb['time']-t))<1e-10]
    result['complete']=bool(not a['failure'] and not b['failure'] and abs(a['completed_end']-a['requested_end'])<1e-10 and abs(b['completed_end']-b['requested_end'])<1e-10)
    result['passed']=False
    if len(pairs)<2:return dict(**result,reason='Insufficient shared recorded times')
    result['shared_start']=float(aa['time'][pairs[0][0]]);result['shared_end']=float(aa['time'][pairs[-1][0]])
    result['errors']={}
    for key in ['gap_arclength_rms','surface_gap_rms']:
        x=np.array([a['diagnostics'][i][key] for i,j in pairs]);y=np.array([b['diagnostics'][j][key] for i,j in pairs])
        result['errors'][key]=dict(maximum_scaled=float(np.max(abs(x-y)/np.maximum(.01,abs(y)))),maximum_absolute=float(np.max(abs(x-y))),final_scaled=float(abs(x[-1]-y[-1])/max(.01,abs(y[-1]))))
    result['energy_passed']=max(a['energy_balance_scaled'],b['energy_balance_scaled'])<1e-3
    result['geometry_passed']=all(r['audit']['minimum_swept_area_ratio']>r['parameters']['area_floor'] and r['audit']['minimum_swept_quality_ratio']>r['parameters']['shape_floor'] and min(r['audit']['minimum_swept_clearance'],r['audit']['minimum_swept_self_clearance'])>0 for r in [a,b])
    result['passed']=bool(result['complete'] and result['refines_both_axes'] and result['energy_passed'] and result['geometry_passed'] and result['shared_end']>=a['requested_end']-1e-10 and all(e['maximum_scaled']<=.05 for e in result['errors'].values()))
    result['input_sha256']=[dict(metadata=digest(p),states=digest(p.with_suffix('.npz'))) for p in [coarse,fine]]
    return result

def select(paths):
    checks=[]
    for coarse,fine in zip(paths,paths[1:]):
        check=compare(coarse,fine);checks.append(check)
        if check['passed']:
            chosen=json.loads(Path(fine).read_text())
            return dict(accepted=True,selected=Path(fine).name,selected_mesh=[chosen['parameters']['angles'],chosen['parameters']['layers']],checks=checks,scope='Agreement of two gap statistics at saved times; not exact error or a fixed-coarse-mesh repair')
    return dict(accepted=False,selected=None,selected_mesh=None,checks=checks,scope='No acceptable resolved result at available levels')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('name');ap.add_argument('--evidence',nargs='+',type=Path);ap.add_argument('--two',action='store_true');ap.add_argument('--end',type=float,default=.2);ap.add_argument('--alpha',type=float,default=16);ap.add_argument('--rho',type=float,default=160);ap.add_argument('--viscosity',type=float,default=.15);ap.add_argument('--seconds',type=float,default=2400);args=ap.parse_args()
    output=ROOT/'data'/(args.name+'-resolution.json')
    if output.exists():raise ValueError('Use a new name; preserve every previous decision')
    if args.evidence:
        result=select(args.evidence)
    else:
        sys.path.insert(0,str(ROOT/'solver'));from experiment import run,Parameters
        paths=[];result=dict(accepted=False,selected=None,checks=[])
        for n,l in [(24,3),(48,5),(72,7)]:
            name=f'{args.name}-{n}x{l}'
            run(name,Parameters(angles=n,layers=l,alpha=args.alpha,rho=args.rho,viscosity=args.viscosity,preload=.35,feedback=args.two),end=args.end,max_seconds=args.seconds)
            paths.append(ROOT/'solver/data'/(name+'.json'))
            if len(paths)>1:result=select(paths)
            if result['accepted']:break
    result['driver_sha256']=digest(__file__)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));raise SystemExit(0 if result['accepted'] else 1)
if __name__=='__main__':main()
