"""Integrity and physical/analysis invariants of the actual completed artifacts."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np
from analyze import geometry,response_metrics
from simulate import dump,tasks

def main():
    p=argparse.ArgumentParser();p.add_argument('--compact',action='store_true');args=p.parse_args();root=Path(__file__).resolve().parent
    rows=list(csv.DictReader((root/'analysis/trajectories.csv').open()));stats=json.loads((root/'analysis/statistics.json').read_text());inv=json.loads((root/'analysis/inventory.json').read_text())
    assert len(rows)==544==len(tasks())==inv['completed'];assert inv['counts']==dict(main=352,long=48,size27=48,wall=48,halfdt=24,friction=24)
    assert stats['main']['mass']['hb_degree']['n']==16
    lookup={r['stem']:r for r in rows};label_errors=[]
    for r in rows:
        for k in ['hb_degree','coordination','nn_exchange','hb_rmst','residence_rmst','temperature']:assert np.isfinite(float(r[k])),(r['stem'],k)
        if r['label']=='True':
            a=lookup[r['stem'].replace('_label','_A')]
            for k in ['hb_degree','coordination','nn_exchange','hb_rmst','residence_rmst','degree_sd']:
                label_errors.append(abs(float(r[k])-float(a[k])))
    assert max(label_errors)<1e-10
    checked=0;max_oh=0.;max_hh=0.;perm_hb=0;response_perm=0.
    if not args.compact:
        expected_hh=2*.09572*np.sin(np.deg2rad(104.52)/2)
        for r in rows:
            path=root/'data'/(r['stem']+'.npz');z=np.load(path);x=z['positions_nm'];assert np.isfinite(x).all();assert len(x)==round(float(r['production_ps'])/.05)+1
            sample=x[::50];oh=np.linalg.norm(sample[:,:,1:]-sample[:,:,0,None,:],axis=-1);hh=np.linalg.norm(sample[:,:,1]-sample[:,:,2],axis=-1)
            max_oh=max(max_oh,float(abs(oh-.09572).max()));max_hh=max(max_hh,float(abs(hh-expected_hh).max()))
            assert abs(z['masses_da'].sum()-(int(r['n'])*18.015324+(2*(2.01410177812-1.007947) if int(r['heavy'])>=0 else 0)))<1e-8
            if r['suite']=='main' and r['label']=='False':
                perm=np.random.default_rng(883+checked).permutation(9);_,h,_=geometry(sample);_,hp,_=geometry(sample[:,perm]);perm_hb+=int(np.any(hp!=h[:,perm][:,:,perm]))
                if 'response_positions_nm' in z:
                    d=dict(z);raw=d['response_positions_nm'];d['response_positions_nm']=raw[:,perm][:,:,:,:,:,perm];d['snapshot_q_nm']=d['snapshot_q_nm'].reshape(3,9,3,3)[:,perm].reshape(3,27,3)
                    v1=response_metrics(z)[0]['response_offdiag'];v2=response_metrics(d)[0]['response_offdiag'];response_perm=max(response_perm,abs(v1-v2))
            checked+=1
        assert max_oh<2e-6 and max_hh<2e-6
        assert perm_hb==0 and response_perm<1e-10
        assert len(list((root/'switch').glob('c[0-9][0-9]_h*.npz')))==160
        assert len(list((root/'refinement').glob('c*.npz')))==72
        manifests={}
        for folder in ['data','switch','refinement']:
            for path in sorted((root/folder).glob('*.npz')):
                manifests[path.relative_to(root).as_posix()]=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        dump(root/'raw_manifest.json',manifests)
    report=(root/'report.html').read_text(encoding='utf-8')
    for text in ['544','0.5 ps','not numerically reliable','not completed','isotope']:
        assert text in report or text in (root/'README.md').read_text(encoding='utf-8'),text
    result=dict(mode='compact' if args.compact else 'full',completed_trajectories=len(rows),raw_checked=checked,raw_files=len(manifests) if not args.compact else None,
                max_group_label_metric_error=max(label_errors),max_OH_constraint_error_nm=max_oh if not args.compact else None,max_HH_constraint_error_nm=max_hh if not args.compact else None,
                heavy_and_light_network_permutation_mismatches=perm_hb if not args.compact else None,max_response_permutation_error=response_perm if not args.compact else None,all_passed=True)
    dump(root/('verification_compact.json' if args.compact else 'verification.json'),result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
