"""Analyze completed source-addition controls, without fitting the q model."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from measure_sources import build, observe, read, save, hashes, original

state=read(HERE/'task-state.json');assert state['status']=='complete' and not state['failures']
ref=read(HERE/'reference-fluid-sources.json')
runs={p.parent.name:read(p) for p in (HERE/'runs').glob('*/result.json')}
assert len(runs)==10
ab,c,d,phi,ops=build()
checks=[]
for name,r in runs.items():
    assert r['status']=='complete' and r['series'][-1]['t']==.5
    assert r['source_hashes']==hashes()
    row=r['series'][-1];u=np.load(HERE/'runs'/name/r['field'])
    observed=observe(u,.5,phi,ops)
    errors={k:abs(observed[k]-row[k]) for k in observed}
    assert max(errors.values())<1e-12
    assert all(np.isfinite(list(x.values())).all() for x in r['series'])
    assert max(x['divergence'] for x in r['series'])<1e-10
    assert r['max_speed_cfl']<.5
    checks.append({'run':name,'endpoint_recalculated_max_absolute_error':max(errors.values()),
       'energy_decreases_at_saved_times':bool(np.all(np.diff([x['energy'] for x in r['series']])<=1e-12)),
       'max_divergence_times_dx':max(x['divergence'] for x in r['series'])*ops[-1],
       'max_high_band_enstrophy_fraction':max(x['high_band_enstrophy_fraction'] for x in r['series'])})
baseline=ref['ab'][-1]['D']
contrasts={}
for case,cn,dn,rn in [('equal','half_c','half_d','ab_cd'),('unequal','seven_c','three_d','ab_uneven')]:
    out={}
    for step in ('base','half'):
        sc=runs[cn+'-'+step]['series'][-1]['S'];sd=runs[dn+'-'+step]['series'][-1]['S']
        mix=ref[rn][-1]['D'] if step=='base' else runs[case+'-half']['series'][-1]['S']
        ec=ref[rn][-1]['E'] if step=='base' else runs[case+'-half']['series'][-1]['E']
        dc,dd=sc-baseline,sd-baseline
        interaction=mix-sc-sd+baseline
        out[step]={'baseline_S':baseline,'C_alone_S':sc,'D_alone_S':sd,'C_contribution':dc,'D_contribution':dd,
           'additive_prediction':baseline+dc+dd,'combined_S':mix,'interaction':interaction,'combined_E':ec,
           'interaction_relative_to_sum_absolute_effects':abs(interaction)/max(abs(dc)+abs(dd),1e-30)}
        assert np.isclose(baseline+dc+dd+interaction,mix,rtol=0,atol=1e-17)
    out['interaction_dt_change']=abs(out['half']['interaction']-out['base']['interaction'])
    out['interaction_dt_change_fraction']=out['interaction_dt_change']/abs(out['half']['interaction']) if abs(out['half']['interaction'])>1e-12 else None
    contrasts[case]=out
mirrorchecks=[]
for step in ('base','half'):
    a=runs['half_c-'+step];b=runs['half_d-'+step]
    ua=np.load(HERE/'runs'/('half_c-'+step)/a['field']);ub=np.load(HERE/'runs'/('half_d-'+step)/b['field'])
    mirror=np.array(original.mirror_of(ua))
    err=float(np.linalg.norm(ub-mirror)/np.linalg.norm(ua))
    assert err<1e-12
    mirrorchecks.append({'step':step,'full_field_relative_mirror_error':err,
       'signed_sum':a['series'][-1]['S']+b['series'][-1]['S'],
       'E_difference':a['series'][-1]['E']-b['series'][-1]['E']})
timechecks=[]
for case in ('half_c','half_d','seven_c','three_d','equal','unequal'):
    rr=runs[case+'-half'];fine=rr['series']
    coarse=ref['ab_cd' if case=='equal' else 'ab_uneven'] if case in ('equal','unequal') else runs[case+'-base']['series']
    assert len(fine)==len(coarse)
    assert all(round(a['t'],3)==round(b['t'],3) for a,b in zip(fine,coarse))
    signedkey='D' if case in ('equal','unequal') else 'S'
    errors={k:max(abs(a[k]-b[signedkey if k=='S' else k]) for a,b in zip(fine,coarse)) for k in ('S','E','W')}
    timechecks.append({'case':case,'coarse_dt':2*rr['dt'],'fine_dt':rr['dt'],'max_absolute_observable_changes':errors})
result={'status':'complete','snapshot_utc':datetime.now(timezone.utc).isoformat(),'n':33,'time':.5,
 'new_runs':10,'reused_source_runs':5,'contrasts':contrasts,'mirror_checks':mirrorchecks,'timestep_checks':timechecks,
 'field_and_record_checks':checks,
 'reading':'Equal source amplitudes cancel the measured signed asymmetry. At70/30, the separate effects predict almost all of the combined endpoint; the interaction contrast is small and stable under dt halving.',
 'limits':['This is an endpoint comparison at33 cubed; high-band enstrophy remains about11.7percent at the endpoint. Spatial convergence is not established.',
 'The contrast measures nonadditivity of source additions relative to AB. It is not a unique fluid-component decomposition or proof of competition between pairings.',
 'The combined runs were already known; this is a retrospective additive comparison, not unseen held-out validation.',
 'No q weights or cubic term were fitted. AB-bias and an independently shaped D were not varied. The full dynamic multi-source model remains uncalibrated.']}
save(HERE/'analysis.json',result)

fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
colors=['#247c8a','#aa6530','#313952']
labels=['C effect','D effect','Combined\nS']
for j,case in enumerate(('equal','unequal')):
    row=contrasts[case]['half']
    vals=[row['C_contribution'],row['D_contribution'],row['combined_S']]
    axs[0,0].bar(np.arange(3)+j*4,vals,color=colors)
axs[0,0].set_xticks([0,1,2,4,5,6],labels+labels)
axs[0,0].tick_params(axis='x',labelsize=9)
axs[0,0].set(title='Signed endpoint: 50/50 (left), 70/30 (right)',ylabel='Signed fluid asymmetry S')
axs[0,0].axhline(0,color='gray',lw=.7)
vals=[contrasts['unequal'][x]['interaction'] for x in ('base','half')]
axs[0,1].bar(['Base step','Half step'],vals,color=['#247c8a','#8eabb2'])
axs[0,1].set(title='70/30 interaction contrast',ylabel='Combined S − additive prediction')
axs[0,1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
for name in ('half_c','half_d','seven_c','three_d','unequal'):
    rows=runs[name+'-half']['series']
    label={'half_c':'0.5 C only','half_d':'0.5 D only','seven_c':'0.7 C only','three_d':'0.3 D only','unequal':'0.7 C + 0.3 D'}[name]
    axs[1,0].plot([x['t'] for x in rows],[x['S'] for x in rows],'.-',label=label)
axs[1,0].set(title='Measured signed asymmetry',xlabel='Model time',ylabel='S');axs[1,0].legend(fontsize=8)
for name in ('equal','unequal'):
    rows=runs[name+'-half']['series'];axs[1,1].plot([x['t'] for x in rows],[x['E'] for x in rows],'.-',label='50/50' if name=='equal' else '70/30')
axs[1,1].set(title='Unsigned mirror difference',xlabel='Model time',ylabel='E');axs[1,1].legend()
for ax in axs.flat:ax.grid(axis='y',alpha=.15)
fig.suptitle('Separate source effects · supplied 33³ fluid · time 0.5',fontsize=15)
fig.savefig(HERE/'source-effects.png',dpi=160);plt.close(fig)
r=contrasts['unequal']['half'];relative=100*abs(r['interaction']/r['combined_S'])
table='| Setting | C contribution | D contribution | Combined S | Interaction contrast |\n| --- | ---: | ---: | ---: | ---: |\n'
for name in ('equal','unequal'):
    v=contrasts[name]['half'];table+='| '+('50/50' if name=='equal' else '70/30')+' | '+' | '.join(f'{v[k]:.10g}' for k in ('C_contribution','D_contribution','combined_S','interaction'))+' |\n'
doc=f'''# Separate C and D contributions

**At 50/50, the measured directions cancel. At 70/30, the separate C and D effects account for all but {relative:.5f}% of the combined signed asymmetry at time 0.5.**

Ten additional fluid runs completed: four missing single-source controls and six half-step checks. Five existing source trajectories were reused. The supplied AB starting field, C, D=mirror(C), 33³ grid, viscosity 0.01, periodic cube side 6, Heun solver and zero force were retained.

![Source contributions and controls](source-effects.png)

## Results

{table}
The table uses half-step endpoints. The 70/30 interaction contrast changes by {100*contrasts['unequal']['interaction_dt_change_fraction']:.5f}% when the timestep is halved. Equal-source runs keep E at roundoff, and the C-only and D-only final velocity fields match after reflection within {max(x['full_field_relative_mirror_error'] for x in mirrorchecks):.3g} relative error.

## What was measured

S is `<u − M[u], phi>`, using the supplied fixed antisymmetric template. It is the signed measurement called D in the original source JSON. **Source D** is C's reflected field. E remains the unsigned normalized mirror difference.

S measures one fixed direction; E measures the full mirror difference relative to the current velocity norm. They can change in different directions over time.

At a fixed endpoint, let S0 be the AB-only result:

```text
C contribution = S(AB + a*C) − S0
D contribution = S(AB + b*D) − S0
interaction = S(AB + a*C + b*D) − S(AB + a*C) − S(AB + b*D) + S0
```

The interaction is the difference from simple addition. These effects depend on the chosen AB background, source weights, observable and endpoint. They do not assign permanent identities to pieces of the evolving fluid.

## What this establishes

For these two source settings, additive effects describe the signed endpoint closely. A large cross-pair term is not needed to explain these measured endpoints. E is not decomposed by adding signed contributions.

This is a retrospective comparison with already recorded combined runs. The full q equation has not been fitted or tested on unseen runs. The cubic term, separate dynamic coupling weights, independent AB bias and independently shaped D remain untested here.

The half-step checks address time integration. The endpoint has about 11.7% of enstrophy near the retained grid cutoff; this experiment does not establish spatial convergence or a continuum mechanism.

## Files and reproduction

[Measured contrasts and verification](analysis.json) · [Protocol](protocol.json) · [Reused source results](reference-fluid-sources.json) · [Run measurements](runs/)

Install the packages in `requirements.txt`, then run `python reproduce.py`. It copies the scripts and reference data into `rerun/`, evolves the controls there, and checks their saved fields. Published measurement records remain available for comparison. The unchanged received source files are in `source/`.
'''
(HERE/'README.md').write_text(doc,encoding='utf-8')
print(json.dumps({'contrasts':contrasts,'mirror_checks':mirrorchecks,'extra_effect_percent_of_combined':relative,'max_cutoff_enstrophy':max(x['max_high_band_enstrophy_fraction'] for x in checks)},indent=2))
