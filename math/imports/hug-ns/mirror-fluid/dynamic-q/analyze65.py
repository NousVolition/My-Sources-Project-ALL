"""Compare new65-grid controls with completed33/49 runs; never refit."""
import hashlib,json
import numpy as np
import analyze as a
import study as s
import matplotlib.pyplot as plt

def main():
    for path,sha in s.read(s.HERE/'training-freeze.json')['files'].items():
        assert hashlib.sha256((s.HERE/path).read_bytes()).hexdigest()==sha,path
    fits=s.read(s.HERE/'fitted-training-models.json');out={}
    fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained')
    for col,(name,*_) in enumerate(s.TEST):
        runs={n:s.read(s.HERE/'runs'/f'{name}-n{n}'/'result.json') for n in (33,49,65)}
        assert all(r['status']=='complete' for r in runs.values())
        r=runs[65];b=s.build(65);u=np.load(s.HERE/'runs'/r['job']['id']/r['field']);last=r['series'][-1]
        checked=s.observe(u,.5,last['q_record'],b)
        errors={key:abs(checked[key]-last[key]) for key in ['q_fluid','E','W','energy','qdot_fluid']}
        assert max(errors.values())<1e-10
        assert r['source_hashes']==s.hashes() and np.isfinite(u).all()
        assert r['max_speed_cfl']<.5
        qold=a.compare(runs[33],runs[49],'q_fluid');qnew=a.compare(runs[49],r,'q_fluid')
        wold=a.compare(runs[33],runs[49],'W');wnew=a.compare(runs[49],r,'W')
        eold=a.compare(runs[33],runs[49],'E');enew=a.compare(runs[49],r,'E')
        energy=np.array([x['energy'] for x in r['series']])
        out[name]=dict(q_33_49=qold,q_49_65=qnew,W_33_49=wold,W_49_65=wnew,E_33_49=eold,E_49_65=enew,
            initial_W={n:v['series'][0]['W'] for n,v in runs.items()},
            max_tail_enstrophy_fraction=max(x['high_band_enstrophy_fraction'] for x in r['series']),
            max_divergence_times_dx=max(x['divergence'] for x in r['series'])*b['ops'][-1],
            max_energy_increase_fraction=float(max(0,np.max(np.diff(energy)))/energy[0]),
            final_field_sha256=hashlib.sha256((s.HERE/'runs'/r['job']['id']/r['field']).read_bytes()).hexdigest(),
            final_recompute_errors=errors,
            predictions={model:{'conditional':a.prediction(r,fit),'future_frozen_inputs':a.prediction(r,fit,.1,True)} for model,fit in fits.items()})
        for n,v in runs.items():
            t=[x['t'] for x in v['series']]
            axes[0,col].plot(t,[x['W'] for x in v['series']],label=f'{n} cubed')
            axes[1,col].plot(t,[x['q_fluid'] for x in v['series']],label=f'{n} cubed')
        axes[0,col].set_title(name.replace('_',' '))
        for ax in axes[:,col]:ax.set_xlabel('Model time');ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Maximum vorticity');axes[1,0].set_ylabel('Signed fluid coordinate');axes[0,0].legend()
    fig.suptitle('Spatial refinement: same four unseen starts',fontsize=17)
    fig.savefig(s.HERE/'grid65-comparison.png',dpi=150);plt.close(fig)
    pooled={model:{mode:a.pooled([v['predictions'][model][mode] for v in out.values()]) for mode in ['conditional','future_frozen_inputs']} for model in a.MODELS}
    result=dict(status='verified',cases=out,scores_on_65=pooled,training_unchanged=True,source_hashes=s.hashes(),
       limitation='No 65-grid timestep repeat. Its base CFL target is the unchanged 0.04; previous half-step controls were at 33. Smaller grid differences do not alone establish spatial convergence.')
    s.save(s.HERE/'grid65-analysis.json',result)
    lines=['# Additional 65-grid check','',
      'Four further runs repeat the same unseen starts at 65 cubed, through time 0.5. The fluid equation, initial formulas, viscosity and frozen training coefficients are unchanged. No 65-grid result is used for fitting.','',
      '| Start | Peak-spin difference 33 to 49 | Peak-spin difference 49 to 65 | Signed-coordinate difference 49 to 65 |','| --- | ---: | ---: | ---: |']
    for name,v in out.items():lines.append(f"| {name} | {100*v['W_33_49']['relative_l2']:.3f}% | {100*v['W_49_65']['relative_l2']:.3f}% | {100*v['q_49_65']['relative_l2']:.5f}% |")
    lines += ['', 'Each value is the relative curve L2 difference against the finer grid.','',
      '![Grid comparison](grid65-comparison.png)','',
      '## Predictions with the original fitted coefficients','',
      '| Model | With measured source histories | With inputs frozen after 0.1 |','| --- | ---: | ---: |']
    for name,v in pooled.items():lines.append(f"| {name.replace('_',' ')} | {100*v['conditional']['relative_l2']:.4f}% | {100*v['future_frozen_inputs']['relative_l2']:.4f}% |")
    lines += ['',result['limitation'],'',
      'All four final fields were independently remeasured. Field hashes, divergence, energy changes and spectral fractions are recorded in [grid65-analysis.json](grid65-analysis.json).', '',
      '[Original 26-run test](README.md) · [Extension protocol](refine65-protocol.json) · [Runner](refine65.py) · [Analyzer](analyze65.py)','']
    (s.HERE/'GRID65.md').write_text('\n'.join(lines),encoding='utf-8')
    finalize_report()
    print(json.dumps({'grid_errors':{k:{'W49_65':v['W_49_65']['relative_l2'],'q49_65':v['q_49_65']['relative_l2']} for k,v in out.items()},'scores65':pooled},indent=2))
def finalize_report():
    p=s.HERE/'grid65-analysis.json';r=s.read(p)
    for case,v in r['cases'].items():
        v['source_features_49_65']=a.compare(s.read(s.HERE/'runs'/f'{case}-n49/result.json'),s.read(s.HERE/'runs'/f'{case}-n65/result.json'),'source_features')
    s.save(p,r)
    old=s.HERE/'GRID65.md';text=old.read_text(encoding='utf-8')
    lead='**Refinement reduces the peak-spin gap from about 28% to 3.16–3.43%. The signed fluid coordinate changes by at most 0.03281%. The fitted full equation still has no consistent advantage over simpler source models, and the cubic term adds no benefit.**\n\n'
    if lead not in text:text=text.replace('# Additional 65-grid check\n\n','# Additional 65-grid check\n\n'+lead,1)
    tail=max(v['max_tail_enstrophy_fraction'] for v in r['cases'].values())
    feat=max(v['source_features_49_65']['relative_l2'] for v in r['cases'].values())
    notes=f'\n## What the refinement changes\n\nThe 65-grid full model gives 1.8908% conditional curve error and 1.5966% frozen-input forecast error. The AB-plus-C model gives slightly smaller errors, 1.8214% and 1.5697%. Model rankings change from the coarse grid; the full equation is not established by these results.\n\nThe maximum source-feature curve change between 49 and 65 is {100*feat:.4f}%. The largest high-band enstrophy fraction on 65 is {100*tail:.5f}%. Energy decreases at saved samples. The starting formulas and their existing x>0 cutoffs were retained; this study does not establish a smooth continuum limit for those formulas.\n\nTo repeat only this extension in a fresh directory, install the existing requirements and run `python reproduce65.py`. The prior 33/49 measurements and frozen fits are used as references; all four 65-grid fields are evolved again.\n'
    if '## What the refinement changes' not in text:text+=notes
    old.write_text(text,encoding='utf-8')
if __name__=='__main__':main()
