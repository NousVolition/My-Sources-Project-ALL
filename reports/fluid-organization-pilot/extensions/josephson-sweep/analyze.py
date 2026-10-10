"""Tables, paired run-level bootstrap intervals and numerical acceptance checks."""
from pathlib import Path
import csv, json
import numpy as np
from model import voltage, from_state
ROOT=Path(__file__).resolve().parent


def canonical(x):
    """Bias-index, mode, alpha, seed, ...; reverse only descending mode."""
    y=x.copy();y[:,2]=x[::-1,2];return y


def save_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    c=json.loads((ROOT/'protocol.json').read_text());out=ROOT/'results';out.mkdir(exist_ok=True)
    n=c['n_independent_initializations'];rng=np.random.default_rng(c['bootstrap_seed'])
    resamples=rng.integers(0,n,(c['bootstrap_resamples'],n))
    def interval(values):
        x=np.asarray(values)
        lo,hi=np.quantile(x[resamples].mean(axis=-1),[.025,.975])
        return {'mean':float(x.mean()),'ci95':[float(lo),float(hi)],'sd':float(x.std(ddof=1))}
    v={};maxpower=0.;raw=[];means=[];valid=True
    for s in c['settings']:
        with np.load(ROOT/f'data/{s["name"]}.npz') as z:
            vi=canonical(z['voltage']);res=canonical(z['residual'])
            valid &= bool(np.isfinite(z['states']).all() and np.all(np.diff(z['dissipated'],axis=1)>=-1e-14))
        maxpower=max(maxpower,float(abs(res).max()));v[s['name']]=vi.sum(-1)
        for j,b in enumerate(c['biases']):
            for m,mode in enumerate(c['modes']):
                for k,a in enumerate(c['alphas']):
                    info=interval(vi[j,m,k].sum(-1))
                    means.append({'setting':s['name'],'dt':s['dt'],'settle':s['settle'],'window':c['base_window'],
                                  'mode':mode,'alpha':a,'bias':b,'mean_total_voltage':info['mean'],
                                  'ci95_low':info['ci95'][0],'ci95_high':info['ci95'][1],'sd':info['sd']})
                    for seed in range(n):
                        raw.append({'setting':s['name'],'mode':mode,'alpha':a,'bias':b,'initialization_id':seed,
                                    'v1':float(vi[j,m,k,seed,0]),'v2':float(vi[j,m,k,seed,1]),
                                    'total_voltage':float(v[s['name']][j,m,k,seed]),
                                    'normalized_power_residual':float(res[j,m,k,seed])})
    windowv={};windowrows=[];blockrows=[]
    with np.load(ROOT/'data/windows.npz') as z:
        p=from_state(z['states'],z['sign']);maxpower=max(maxpower,float(abs(z['residual']).max()))
        valid &= bool(np.isfinite(p).all() and np.all(np.diff(z['dissipated'],axis=0)>=-1e-14))
        for window in c['window_checks']:
            vi=canonical(voltage(p[0],p[round(window/c['record_stride'])],window));windowv[int(window)]=vi.sum(-1)
            for j,b in enumerate(c['biases']):
                for m,mode in enumerate(c['modes']):
                    for k,a in enumerate(c['alphas']):
                        for seed in range(n):
                            windowrows.append({'window':window,'mode':mode,'alpha':a,'bias':b,'initialization_id':seed,
                                               'v1':float(vi[j,m,k,seed,0]),'v2':float(vi[j,m,k,seed,1]),
                                               'total_voltage':float(vi[j,m,k,seed].sum())})
        for block in range(4):
            every=round(200/c['record_stride']);vi=canonical(voltage(p[block*every],p[(block+1)*every],200))
            for j,b in enumerate(c['biases']):
                for m,mode in enumerate(c['modes']):
                    for k,a in enumerate(c['alphas']):
                        for seed in range(n):
                            blockrows.append({'block_start':block*200,'mode':mode,'alpha':a,'bias':b,'initialization_id':seed,
                                              'total_voltage':float(vi[j,m,k,seed].sum())})
    replayerror=float(abs(windowv[200]-v['both']).max())
    critical=(np.array(c['biases'])>=c['critical_window'][0])&(np.array(c['biases'])<=c['critical_window'][1])
    step1=abs(v['half_step']-v['base']);step2=abs(v['both']-v['double_settle'])
    settling=abs(v['both']-v['half_step']);windowdelta=abs(windowv[800]-windowv[200])
    regimes=[]
    for k,a in enumerate(c['alphas']):
        gap200=abs(windowv[200][critical,1,k]-windowv[200][critical,2,k]).mean(axis=0)
        gap800=abs(windowv[800][critical,1,k]-windowv[800][critical,2,k]).mean(axis=0)
        regimes.append({'alpha':a,
            'max_step_change':float(max(step1[:,:,k].max(),step2[:,:,k].max())),
            'max_settling_change':float(settling[:,:,k].max()),
            'max_window_change':float(windowdelta[:,:,k].max()),
            'critical_mean_abs_up_down_gap_W200':interval(gap200),
            'critical_mean_abs_up_down_gap_W800':interval(gap800),
            'paired_gap_change_W800_minus_W200':interval(gap800-gap200),
            'max_abs_subcritical_voltage_both':float(abs(v['both'][np.array(c['biases'])<1,:,k]).max())})
    ref=json.loads((out/'independent_validation.json').read_text())
    metrics={'max_absolute_voltage_step_difference':float(max(step1.max(),step2.max())),
             'max_normalized_power_balance_residual':maxpower,
             'max_reference_voltage_error':ref['max_reference_voltage_error'],
             'max_analytic_cycle_relative_error':ref['max_analytic_cycle_relative_error']}
    gates={k:{'value':val,'limit':c['gates'][k],'passed':val<c['gates'][k]} for k,val in metrics.items()}
    gates['finite_and_nonnegative_dissipation']={'passed':valid}
    gates['window_replay_reproduces_W200']={'value':replayerror,'limit':1e-11,'passed':replayerror<1e-11}
    counts={'primary_segments':len(raw),'independent_initial_draws':n,'voltage_window_estimates':len(windowrows),
            'long_window_replays':len(windowrows)//3,'independent_solver_runs':len(ref['independent_solver_runs']),
            'analytic_cycle_runs':len(ref['analytic_cycles'])}
    summary={'counts':counts,'regimes':regimes,'numerical_gates':gates,'all_gates_passed':all(g['passed'] for g in gates.values()),
             'critical_definition':'Average absolute up-minus-down total voltage across sampled currents 0.98 to 1.02, computed per initial phase draw before bootstrap.',
             'interval_scope':'Paired percentile bootstrap over 12 independent initial phase draws, preserving all currents/modes/settings within a draw; not a population confidence interval for real devices.'}
    for name,rows in [('primary_voltages',raw),('primary_means',means),('window_voltages',windowrows),('disjoint_blocks',blockrows)]:
        save_csv(out/f'{name}.csv',rows)
    flat=[]
    for r in regimes:
        row={}
        for k,val in r.items():
            if isinstance(val,dict):
                row[k]=val['mean'];row[k+'_ci95_low']=val['ci95'][0];row[k+'_ci95_high']=val['ci95'][1]
            else:row[k]=val
        flat.append(row)
    save_csv(out/'regime_summary.csv',flat)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'verification.json').write_text(json.dumps({'all_passed':summary['all_gates_passed'],'gates':gates},indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not summary['all_gates_passed']:
        print('Original fixed-step gate failed. Retain this finding; run refine.py and analyze_reference.py before interpretation.')


if __name__=='__main__':main()
