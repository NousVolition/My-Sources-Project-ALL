"""Summarize resolved full-history references separately from original checks."""
from pathlib import Path
import json
import numpy as np
from analyze import canonical,save_csv
from model import from_state,voltage
ROOT=Path(__file__).resolve().parent


def main():
    c=json.loads((ROOT/'protocol.json').read_text());refinement=json.loads((ROOT/'results/refinement.json').read_text())
    if not refinement['all_passed']:raise RuntimeError('Refinement failed')
    n=c['n_independent_initializations'];rng=np.random.default_rng(c['bootstrap_seed'])
    rs=rng.integers(0,n,(c['bootstrap_resamples'],n));raw=[];means=[];reference={};maxpower=0
    def ci(x):
        lo,hi=np.quantile(x[rs].mean(-1),[.025,.975])
        return {'mean':float(x.mean()),'ci95':[float(lo),float(hi)]}
    for settle in [100,200]:
        with np.load(ROOT/f'data/reference_S{settle}.npz') as z:
            vi=canonical(z['voltage']);maxpower=max(maxpower,float(abs(z['residual']).max()))
        reference[settle]=vi.sum(-1)
        for j,b in enumerate(c['biases']):
            for m,mode in enumerate(c['modes']):
                for k,a in enumerate(c['alphas']):
                    stats=ci(vi[j,m,k].sum(-1))
                    means.append({'settle':settle,'window':200,'mode':mode,'alpha':a,'bias':b,
                                  'mean_total_voltage':stats['mean'],'ci95_low':stats['ci95'][0],'ci95_high':stats['ci95'][1]})
                    for seed in range(n):raw.append({'settle':settle,'window':200,'mode':mode,'alpha':a,'bias':b,'initialization_id':seed,
                        'v1':float(vi[j,m,k,seed,0]),'v2':float(vi[j,m,k,seed,1]),'total_voltage':float(vi[j,m,k,seed].sum())})
    wv={};windowrows=[];blockrows=[]
    with np.load(ROOT/'data/reference_windows.npz') as z:
        p=from_state(z['states'],z['sign']);maxpower=max(maxpower,float(abs(z['residual']).max()))
    for w in [200,400,800]:
        vi=canonical(voltage(p[0],p[w//5],w));wv[w]=vi.sum(-1)
        for j,b in enumerate(c['biases']):
            for m,mode in enumerate(c['modes']):
                for k,a in enumerate(c['alphas']):
                    for seed in range(n):windowrows.append({'window':w,'mode':mode,'alpha':a,'bias':b,'initialization_id':seed,
                        'v1':float(vi[j,m,k,seed,0]),'v2':float(vi[j,m,k,seed,1]),'total_voltage':float(vi[j,m,k,seed].sum())})
    for block in range(4):
        vi=canonical(voltage(p[block*40],p[(block+1)*40],200))
        for j,b in enumerate(c['biases']):
            for m,mode in enumerate(c['modes']):
                for k,a in enumerate(c['alphas']):
                    for seed in range(n):blockrows.append({'block_start':block*200,'mode':mode,'alpha':a,'bias':b,'initialization_id':seed,
                                                         'total_voltage':float(vi[j,m,k,seed].sum())})
    mask=(np.array(c['biases'])>=.98)&(np.array(c['biases'])<=1.02)
    with np.load(ROOT/'data/both.npz') as z:original=canonical(z['voltage']).sum(-1)
    regimes=[]
    for k,a in enumerate(c['alphas']):
        gap200=abs(wv[200][mask,1,k]-wv[200][mask,2,k]).mean(0)
        gap800=abs(wv[800][mask,1,k]-wv[800][mask,2,k]).mean(0)
        regimes.append({'alpha':a,'original_both_vs_reference_max_error':float(abs(original[:,:,k]-reference[200][:,:,k]).max()),
                        'max_settling_change':float(abs(reference[200][:,:,k]-reference[100][:,:,k]).max()),
                        'max_window_change':float(abs(wv[800][:,:,k]-wv[200][:,:,k]).max()),
                        'critical_mean_abs_up_down_gap_W200':ci(gap200),'critical_mean_abs_up_down_gap_W800':ci(gap800),
                        'paired_gap_change_W800_minus_W200':ci(gap800-gap200)})
    result={'regimes':regimes,'max_normalized_power_residual':maxpower,
            'power_gate_passed':maxpower<c['gates']['max_normalized_power_balance_residual'],
            'reference_window_replay_W200_max_difference':float(abs(reference[200]-wv[200]).max()),
            'refinement':refinement,'primary_reference_segments':len(raw),'targeted_refinement_segments':612,
            'interpretation':'Original fixed-step comparisons remain recorded, including their failed gate. These separately labeled reference summaries use full-history DOP853 replays; no averaging window is mixed with another sweep history.'}
    flat=[]
    for r in regimes:
        row={}
        for k,v in r.items():
            if isinstance(v,dict):row.update({k:v['mean'],k+'_ci95_low':v['ci95'][0],k+'_ci95_high':v['ci95'][1]})
            else:row[k]=v
        flat.append(row)
    for name,rows in [('reference_voltages',raw),('reference_means',means),('reference_windows',windowrows),('reference_disjoint_blocks',blockrows),('reference_regimes',flat)]:
        save_csv(ROOT/f'results/{name}.csv',rows)
    (ROOT/'results/reference_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['power_gate_passed'] or result['reference_window_replay_W200_max_difference']>2e-6:raise RuntimeError('Reference consistency failed')


if __name__=='__main__':main()
