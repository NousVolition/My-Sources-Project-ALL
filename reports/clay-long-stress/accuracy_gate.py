"""Fail closed when a mesh or integration comparison has not met its target.

Usage: python accuracy_gate.py coarse.json fine.json [--kind mesh|time]
Neither the finest run nor two agreeing coarse runs are assumed exact.
"""
import argparse,json
from pathlib import Path
import numpy as np
def assess(coarse,fine,kind='mesh'):
    a=json.loads(Path(coarse).read_text());b=json.loads(Path(fine).read_text())
    result=dict(kind=kind,coarse=str(coarse),fine=str(fine),complete=not a['failure'] and not b['failure'])
    pa=dict(a['parameters']);pb=dict(b['parameters'])
    if kind=='mesh':
        for k in ['angles','layers']:pa.pop(k);pb.pop(k)
        result['refines_mesh']=b['parameters']['angles']>a['parameters']['angles'] and b['parameters']['layers']>a['parameters']['layers']
    assert pa==pb,'Comparisons must use the same physical model and parameters'
    assert a['initial']==b['initial'] and a['requested_end']==b['requested_end'],'Initial state and requested interval must match'
    aa=np.load(Path(coarse).with_suffix('.npz'));bb=np.load(Path(fine).with_suffix('.npz'))
    pairs=[(i,int(np.argmin(abs(bb['time']-t)))) for i,t in enumerate(aa['time']) if np.min(abs(bb['time']-t))<1e-10]
    if len(pairs)<2:return dict(**result,passed=False,reason='Insufficient shared saved times')
    ia,ib=np.array(pairs).T
    result['shared_end']=float(aa['time'][ia[-1]])
    if kind=='mesh':
        ga=np.array([a['diagnostics'][i]['gap_arclength_rms'] for i in ia]);gb=np.array([b['diagnostics'][i]['gap_arclength_rms'] for i in ib])
        difference=abs(ga-gb);scaled=difference/np.maximum(.01,abs(gb))
        result.update(final_gap_relative=float(scaled[-1]),maximum_gap_scaled_difference=float(np.max(scaled)),maximum_gap_absolute_difference=float(np.max(difference)),gap_scale_floor=.01,tolerance=.05)
        agreed=bool(np.max(scaled)<=.05 and result['refines_mesh'])
    else:
        ya=aa['state'][ia,:-4];yb=bb['state'][ib,:-4];error=float(np.max(abs(ya-yb)/np.maximum(1,np.max(abs(yb),axis=0))))
        result.update(scaled_state_error=error,tolerance=1e-4);agreed=error<=1e-4
    result['passed']=bool(result['complete'] and result['shared_end']>=a['requested_end']-1e-10 and agreed)
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('coarse');ap.add_argument('fine');ap.add_argument('--kind',choices=['mesh','time'],default='mesh');ap.add_argument('--output');args=ap.parse_args()
    result=assess(args.coarse,args.fine,args.kind);text=json.dumps(result,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text);raise SystemExit(0 if result['passed'] else 1)
