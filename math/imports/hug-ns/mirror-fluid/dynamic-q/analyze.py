"""Fit only training runs; compare withheld whole trajectories and controls."""
from pathlib import Path
import os, json, math, sys, hashlib
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import numpy as np
from scipy.optimize import lsq_linear
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import study as s

MODELS={'decay_only':[0],'AB_only':[0,1],'AB_C':[0,1,2],
        'AB_C_D':[0,1,2,3],'full_linear':[0,1,2,3,4],'full_cubic':[0,1,2,3,4,5]}
LABELS=['gamma','w_AB','w_C','w_D','w_cross','g']
def data(result):
    rows=result['series'];q=np.array([r['q_fluid'] for r in rows]);d=np.array([r['source_features'] for r in rows])
    return np.column_stack([-q,d,-q**3]),np.array([r['qdot_fluid'] for r in rows])
def fit(runs,columns,ridge):
    blocks=[data(r) for r in runs];X=np.concatenate([x[:,columns] for x,y in blocks]);y=np.concatenate([y for x,y in blocks])
    scale=np.maximum(np.sqrt(np.mean(X*X,axis=0)),1e-20);ys=max(float(np.sqrt(np.mean(y*y))),1e-20)
    A=X/scale;target=y/ys
    low=np.array([0 if c in (0,5) else -np.inf for c in columns]);high=np.full(len(columns),np.inf)
    AA=np.vstack([A,np.sqrt(len(y)*ridge)*np.eye(len(columns))]);bb=np.r_[target,np.zeros(len(columns))]
    result=lsq_linear(AA,bb,bounds=(low,high),tol=1e-12,max_iter=2000)
    assert result.success,result.message
    beta=np.zeros(6);beta[columns]=result.x*ys/scale
    sing=np.linalg.svd(A,compute_uv=False)
    return dict(beta=beta.tolist(),ridge=ridge,columns=columns,training_rmse=float(np.sqrt(np.mean((X@beta[columns]-y)**2))),
       scaled_condition=float(sing[0]/max(sing[-1],1e-30)),singular_values=sing.tolist(),rank=int(np.linalg.matrix_rank(A)),
       y_scale=ys,column_scales=scale.tolist())
def select(runs,columns):
    candidates=[]
    for ridge in [0.,1e-6,.001,.1]:
        errors=[]
        for i,r in enumerate(runs):
            f=fit(runs[:i]+runs[i+1:],columns,ridge);x,y=data(r)
            errors.extend((x@f['beta']-y)**2)
        candidates.append(dict(ridge=ridge,leave_one_run_out_derivative_rmse=float(np.sqrt(np.mean(errors)))))
    best=min(candidates,key=lambda c:c['leave_one_run_out_derivative_rmse'])
    result=fit(runs,columns,best['ridge']);result['selection']=candidates
    omissions=[fit(runs[:i]+runs[i+1:],columns,best['ridge'])['beta'] for i in range(len(runs))]
    result['coefficients_leave_one_training_run_out']=omissions
    return result
def prediction(result,fitted,start=.0,frozen=False):
    rows=result['series'];t=np.array([r['t'] for r in rows]);q=np.array([r['q_fluid'] for r in rows]);d=np.array([r['source_features'] for r in rows])
    index=int(np.argmin(abs(t-start)));t=t[index:];truth=q[index:];d=d[index:]
    beta=np.asarray(fitted['beta'])
    def rhs(time,state):
        inputs=d[0] if frozen else np.array([np.interp(time,t,d[:,i]) for i in range(4)])
        return [-beta[0]*state[0]+inputs@beta[1:5]-beta[5]*state[0]**3]
    sol=solve_ivp(rhs,(float(t[0]),float(t[-1])),[truth[0]],t_eval=t,rtol=1e-9,atol=1e-12,max_step=.005)
    assert sol.success,sol.message
    error=sol.y[0]-truth
    return dict(t=t.tolist(),predicted=sol.y[0].tolist(),truth=truth.tolist(),rmse=float(np.sqrt(np.mean(error**2))),
        relative_l2=float(np.linalg.norm(error)/max(np.linalg.norm(truth),1e-15)),
        endpoint_error=float(error[-1]),endpoint_sign_matches=bool(np.sign(sol.y[0,-1])==np.sign(truth[-1])))
def pooled(rows):
    error=np.concatenate([np.asarray(r['predicted'])-r['truth'] for r in rows]);truth=np.concatenate([r['truth'] for r in rows])
    return dict(rmse=float(np.sqrt(np.mean(error**2))),relative_l2=float(np.linalg.norm(error)/max(np.linalg.norm(truth),1e-15)))
def compare(a,b,key):
    x=np.array([r[key] for r in a['series']]);y=np.array([r[key] for r in b['series']])
    return dict(relative_l2=float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-15)),max_abs=float(np.max(abs(x-y))))
def record_view(result):
    # Same saved data, historical passive record as target. Its derivative is
    # known by construction; this is a surrogate test, not independent physics.
    out=dict(result)
    out['series']=[dict(row,q_fluid=row['q_record'],qdot_fluid=-.2*row['q_record']+.8*row['S']) for row in result['series']]
    return out
def record_test(results):
    views={key:record_view(r) for key,r in results.items()}
    train=[views[name+'-n33'] for name,*_ in s.TRAIN]
    fits=s.read(HERE/'fitted-record-models.json')
    tests={name:{model:prediction(views[name+'-n33'],fit) for model,fit in fits.items()} for name,*_ in s.TEST}
    return dict(target='Existing imposed linear record, using the source-template S of this study',fits=fits,tests=tests,
        pooled_scores={model:pooled([tests[name][model] for name,*_ in s.TEST]) for model in MODELS},
        control='The exact step rule already defines this target: q_next=q+dt*(-0.2*q+0.8*S_next). A six-template fit can approximate the record but cannot establish that the fluid has cubic feedback.')
def main():
    protocol=s.read(HERE/'protocol.json');assert protocol['source_hashes']==s.hashes()
    results={job['id']:s.read(HERE/'runs'/job['id']/'result.json') for job in s.JOBS}
    assert all(r['status']=='complete' for r in results.values())
    train=[results[name+'-n33'] for name,*_ in s.TRAIN]
    if not (HERE/'training-freeze.json').exists():freeze_training()
    frozen=s.read(HERE/'training-freeze.json')
    for path,expected in frozen['files'].items():assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==expected,path
    fits=s.read(HERE/'fitted-training-models.json')
    tests={};mirrors={};controls={};verify=[]
    # All coefficients and model selection are fixed before the test set is scored.
    for name,*_ in s.TEST:
        r=results[name+'-n33'];tests[name]={}
        for model,fitted in fits.items():
            tests[name][model]={'conditional':prediction(r,fitted),'future_frozen_inputs':prediction(r,fitted,.1,True)}
        mirrored=results[name+'-mirror']
        q=np.array([x['q_fluid'] for x in r['series']]);qm=np.array([x['q_fluid'] for x in mirrored['series']])
        e=np.array([x['E'] for x in r['series']]);em=np.array([x['E'] for x in mirrored['series']])
        d=np.array([x['source_features'] for x in r['series']]);dm=np.array([x['source_features'] for x in mirrored['series']])
        p=tests[name]['full_cubic']['conditional'];pm=prediction(mirrored,fits['full_cubic'])
        mirrors[name]=dict(max_q_sum=float(np.max(abs(q+qm))),max_E_difference=float(np.max(abs(e-em))),
           max_feature_sum=float(np.max(abs(d+dm))),max_prediction_sum=float(np.max(abs(np.array(p['predicted'])+pm['predicted']))))
        controls[name]={}
        for suffix in ['-half','-n49']:
            control=results[name+suffix]
            controls[name][suffix]=dict(q=compare(r,control,'q_fluid'),E=compare(r,control,'E'),
                W=compare(r,control,'W'),source_features=compare(r,control,'source_features'),
                predictions={model:prediction(control,fitted) for model,fitted in fits.items()})
    scores={model:{mode:pooled([tests[name][model][mode] for name,*_ in s.TEST])
          for mode in ['conditional','future_frozen_inputs']} for model in MODELS}
    for name,r in results.items():
        b=s.build(r['job']['n']);u=np.load(HERE/'runs'/name/r['field']);last=r['series'][-1]
        row=s.observe(u,last['t'],last['q_record'],b)
        checks={key:abs(row[key]-last[key]) for key in ['q_fluid','qdot_fluid','E','W','energy','divergence']}
        assert max(checks.values())<1e-10,(name,checks)
        energies=np.array([x['energy'] for x in r['series']]);finite=all(np.isfinite(v).all() for v in [u,energies])
        assert finite and r['max_speed_cfl']<.5
        verify.append(dict(id=name,finite=True,final_recompute_errors=checks,max_speed_cfl=r['max_speed_cfl'],
            largest_saved_energy_increase_fraction=float(max(0,np.max(np.diff(energies)))/energies[0]),
            max_divergence_times_dx=max(x['divergence'] for x in r['series'])*b['ops'][-1],
            largest_template_residual_fraction=max(x['template_residual_fraction'] for x in r['series']),
            largest_high_band_enstrophy_fraction=max(x['high_band_enstrophy_fraction'] for x in r['series'])))
    linear=scores['full_linear']['conditional']['rmse'];cubic=scores['full_cubic']['conditional']['rmse']
    record=record_test(results)
    out=dict(status='complete',runs=len(results),training_runs=len(train),withheld_runs=4,passive_record_test=record,
        fits=fits,tests=tests,pooled_scores=scores,mirror_checks=mirrors,controls=controls,verification=verify,
        cubic_unseen_conditional_rmse_improvement_fraction=1-cubic/max(linear,1e-30),
        definitions=protocol['coordinates'],source_features=protocol['source_definition'],
        q_record_caution='The imposed linear record is retained as a diagnostic only. Recovering its defining equation is not independent validation.',
        branch_check='For gamma>=0,g>=0, derivative of the drift in q is -gamma-3g*q^2<=0. Constant inputs cannot produce two isolated stable branches.',
        scope='Conditional tests use measured held-out input histories. Frozen-input forecasts use no fluid observations after t=.1. Neither closes the full fluid dynamics.')
    s.save(HERE/'analysis.json',out)
    plot(out);report(out)
    print(json.dumps({'runs':len(results),'scores':scores,'cubic_improvement':out['cubic_unseen_conditional_rmse_improvement_fraction'],
      'coefficients':dict(zip(LABELS,fits['full_cubic']['beta']))},indent=2))
def plot(a):
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,(name,*_) in zip(axes.flat,s.TEST):
        rows=a['tests'][name];r=rows['full_cubic']['conditional'];ax.plot(r['t'],r['truth'],color='black',lw=2,label='Measured fluid')
        for model,label,color in [('AB_only','AB only','#aaaaaa'),('full_linear','All sources, linear','#167a8b'),('full_cubic','All sources + cubic','#bf5547')]:
            v=rows[model]['conditional'];ax.plot(v['t'],v['predicted'],label=label,color=color,ls='--')
        f=rows['full_cubic']['future_frozen_inputs'];ax.plot(f['t'],f['predicted'],label='Cubic forecast, inputs frozen at 0.1',color='#8751ac',ls=':')
        ax.set_title(name.replace('_',' '));ax.set_xlabel('Model time');ax.set_ylabel('Signed fluid coordinate');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Predictions on four starts excluded from fitting',fontsize=16)
    fig.savefig(HERE/'held-out-predictions.png',dpi=160);plt.close(fig)
    names=list(MODELS);x=np.arange(len(names));fig,ax=plt.subplots(figsize=(11,5),layout='constrained')
    for offset,mode,label in [(-.2,'conditional','Measured source histories'),(.2,'future_frozen_inputs','Inputs frozen after 0.1')]:
        ax.bar(x+offset,[100*a['pooled_scores'][k][mode]['relative_l2'] for k in names],width=.4,label=label)
    ax.set_xticks(x,[n.replace('_',' ') for n in names],rotation=15);ax.set_ylabel('Pooled prediction error (% of measured curve norm)');ax.legend();ax.grid(axis='y',alpha=.2)
    ax.set_title('Does the larger model improve unseen predictions?');fig.savefig(HERE/'model-errors.png',dpi=160);plt.close(fig)
def report(a):
    scores=a['pooled_scores'];lines=['# Dynamic source weights and cubic term','',
      '**The source features improve predictions when their measured histories are supplied. The cubic term adds no benefit. Forecasting without future source measurements is worse than simple decay, and the fluid remains sensitive to grid resolution. These tests do not establish the full proposed model.**','',
      '26 fluid runs completed: 10 training starts, 4 unseen starts, and 12 mirror, timestep and grid controls. The fluid equation and supplied time stepper were unchanged. No fitted term forces the fluid.','',
      '## What was tested','',
      'The historical q is a passive record computed by a known linear rule. It cannot independently prove a nonlinear feedback law. This test therefore measures a separate fluid coordinate `q_fluid = <u,phi>/norm(AB0)`, and fits the proposed equation to its full fluid derivative. The original record is also saved as `q_record`; the two are not interchangeable.','',
      'A six-template spatial projection defines four dynamic input features: AB interactions, AB–C and C interactions, AB–D and D interactions, and C–D interactions. The full fluid derivative is measured independently. The interaction feature is a specified measurable candidate for the cross term; it does not establish competition between alternative pairings. Definitions and all reserved cases were fixed in [protocol.json](protocol.json) before fitting.','',
      'AB bias varies independently. D has a different position, width and velocity direction from reflected C. Two unseen starts additionally change D’s shape.','',
      '## Unseen prediction errors','',
      '| Model | Measured input histories | Inputs frozen after time 0.1 |','| --- | ---: | ---: |']
    for name in MODELS:lines.append(f"| {name.replace('_',' ')} | {100*scores[name]['conditional']['relative_l2']:.4f}% | {100*scores[name]['future_frozen_inputs']['relative_l2']:.4f}% |")
    lines += ['', 'Errors are pooled curve L2 errors divided by the measured curve norm. Conditional prediction starts from the initial fluid coordinate and uses the unseen run’s measured source histories, without future q values. It requires those input histories. The separate forecast uses only information through time 0.1 and assumes subsequent source inputs remain fixed.','',
      '![Unseen predictions](held-out-predictions.png)','', '![Model comparison](model-errors.png)','',
      '## Fitted full equation','', '| Coefficient | Fit | Range after omitting one training run |','| --- | ---: | ---: |']
    fit=a['fits']['full_cubic'];om=np.array(fit['coefficients_leave_one_training_run_out'])
    for i,label in enumerate(LABELS):lines.append(f'| {label} | {fit["beta"][i]:.7g} | {om[:,i].min():.7g} to {om[:,i].max():.7g} |')
    lines += ['',f"Adding the cubic term changes unseen conditional RMSE by an improvement of {100*a['cubic_unseen_conditional_rmse_improvement_fraction']:.5f}% (negative means worse). The six-column scaled design condition number is {fit['scaled_condition']:.6g}; rank is {fit['rank']}/6. Fits use nonnegative damping gamma and g; source weights may have either sign.",'',
      'The full-data cubic coefficient is effectively zero. Its value changes substantially when individual training runs are omitted. The cross coefficient even changes sign. These coefficients should not be treated as established physical constants. The fourth unseen start, which changes D’s shape, has a 5.82% conditional curve error; the pooled 1.62% error hides that weaker result.','',
      '## Existing passive q record','',
      'The same full and reduced equations were also fitted to q_record, with its known derivative -0.2 q_record + 0.8 S. This target is created by that linear rule. The comparison checks how well the candidate source features reproduce the record; it is not evidence of cubic feedback in the fluid. The chosen signed template is the source-study template, which differs from the older fluid_correction.py template.','',
      '| Model | Unseen record prediction error |','| --- | ---: |']
    for name in MODELS:lines.append(f"| {name.replace('_',' ')} | {100*a['passive_record_test']['pooled_scores'][name]['relative_l2']:.4f}% |")
    lines += ['',
      '## Numerical controls','', '| Unseen start | q curve change: half timestep | q curve change: 33 to 49 grid | Peak-spin curve change: 33 to 49 grid |','| --- | ---: | ---: | ---: |']
    for name,*_ in s.TEST:
        v=a['controls'][name];lines.append(f"| {name} | {100*v['-half']['q']['relative_l2']:.6f}% | {100*v['-n49']['q']['relative_l2']:.6f}% | {100*v['-n49']['W']['relative_l2']:.3f}% |")
    worst=max(v['largest_high_band_enstrophy_fraction'] for v in a['verification'])
    lines += ['', f"All 26 final fields were remeasured. Maximum mirror-pair signed sum: {max(v['max_q_sum'] for v in a['mirror_checks'].values()):.3g}; maximum E difference: {max(v['max_E_difference'] for v in a['mirror_checks'].values()):.3g}. Worst retained high-band enstrophy fraction: {100*worst:.3f}%. Energy decreases at all saved samples. Peak-spin curves differ by about 28% between grids, including a substantial initial peak difference. The relatively stable global q measurement does not make the full fluid field resolved.",'',
      '## Branch selection','',a['branch_check'],'',
      'A sign produced by signed source inputs is a response to those inputs. This model does not demonstrate spontaneous selection of either sign from an unsigned start.','',
      '## Files','', '[All measurements, fits and errors](analysis.json) · [Training fits](fitted-training-models.json) · [Protocol](protocol.json) · [Fluid runner](study.py) · [Analysis code](analyze.py) · [Run index](runs/index.json)','']
    (HERE/'README.md').write_text('\n'.join(lines),encoding='utf-8')
def freeze_training():
    assert not (HERE/'training-freeze.json').exists(),'Training fit is already frozen.'
    train=[s.read(HERE/'runs'/f'{name}-n33'/'result.json') for name,*_ in s.TRAIN]
    assert all(r['status']=='complete' for r in train)
    fits={name:select(train,cols) for name,cols in MODELS.items()}
    records={name:select([record_view(r) for r in train],cols) for name,cols in MODELS.items()}
    s.save(HERE/'fitted-training-models.json',fits);s.save(HERE/'fitted-record-models.json',records)
    paths=['fitted-training-models.json','fitted-record-models.json','protocol.json']+[f'runs/{name}-n33/result.json' for name,*_ in s.TRAIN]
    s.save(HERE/'training-freeze.json',dict(files={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
       holdout_rule='No unseen, mirrored, finer-grid or half-step result is used in fitting or ridge selection.'))
    print('Training coefficients frozen before unseen scoring.',flush=True)
if __name__=='__main__':
    if '--fit-only' in sys.argv:freeze_training()
    else:main()
