"""Create figures and derived checks from saved raw outputs; never reruns physics."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import bootstrap_mean,hazard_model
from run_study import write_csv,write_json

ROOT=Path(__file__).resolve().parent


def read(name):
    with open(ROOT/'data'/name,encoding='utf-8') as f: return list(csv.DictReader(f))


def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=160,bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)


def main():
    (ROOT/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.facecolor':'white','axes.grid':True,'grid.alpha':.2})
    summary=json.loads((ROOT/'summary.json').read_text())
    r=read('solute_timeseries.csv'); p=read('particle_timeseries.csv')
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for flow,color in [('still','#718096'),('stirred','#167d9a')]:
        vals=[]
        for seed in range(100,108):
            q=[a for a in r if a['flow']==flow and int(a['seed'])==seed]
            t=np.array([float(a['t']) for a in q]); vals.append([float(a['variance'])/.10125 for a in q])
        vals=np.array(vals)
        ax[0].plot(t,vals.mean(axis=0),color=color,label=flow)
        ax[0].fill_between(t,vals.min(axis=0),vals.max(axis=0),color=color,alpha=.15)
    ax[0].axhline(.05,c='black',ls=':',label='95% variance reduction')
    ax[0].set(xlabel='Time (s, illustrative mapping)',ylabel='Solute variance / initial variance',yscale='log',title='Dissolved solute: 8 flow seeds')
    ax[0].legend(fontsize=8)
    for mode,color in [('inert','#167d9a'),('brownian','#cc6e34'),('motile_assumed','#7449a3')]:
        vals=[]
        for seed in range(100,108):
            q=[a for a in p if a['mode']==mode and int(a['seed'])==seed]
            t=[float(a['t']) for a in q]; vals.append([float(a['cv2_minus_uniform_noise']) for a in q])
        ax[1].plot(t,np.mean(vals,axis=0),label=mode.replace('_',' '),color=color)
    ax[1].set(xlabel='Time (s)',ylabel='Coarse particle CV² minus noise floor',title='Particles: fixed 16 × 16 observation bins')
    ax[1].legend(fontsize=8)
    fig.suptitle('COMPUTATIONAL PILOT · mixing depends on scale and transport',fontsize=13)
    save(fig,'mixing')

    fig,ax=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    m=np.load(ROOT/'data/solute_stirred_maps.npz')
    q=np.load(ROOT/'data/particles_inert_100.npz')
    for j,i in enumerate([0,20,80]):
        im=ax[0,j].imshow(m['concentration'][i].T,origin='lower',extent=(0,1,0,1),vmin=0,vmax=1,cmap='viridis')
        ax[0,j].set(title=f'Solute, t = {m["time"][i]:g} s',xlabel='x / L',ylabel='z / L')
        ip=ax[1,j].imshow((q['counts'][i]*256/8192).T,origin='lower',extent=(0,1,0,1),vmin=0,vmax=2,cmap='magma')
        ax[1,j].set(title=f'Inert particles, t = {q["time"][i]:g} s',xlabel='x / L',ylabel='z / L')
        ax[0,j].grid(False); ax[1,j].grid(False)
    fig.colorbar(im,ax=ax[0,:],label='Source-water fraction (dimensionless)',shrink=.8)
    fig.colorbar(ip,ax=ax[1,:],label='Binned particle density / mean',shrink=.8)
    fig.suptitle('SIMULATED concentration maps · seed 100 · no particle feedback',fontsize=13)
    save(fig,'maps')

    sr=read('stratification_timeseries.csv')
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for j,mech in enumerate(['density','thermal']):
        for N2,label,color in [(0,'neutral','#718096'),(4,'stable','#167d9a'),(-4,'unstable','#cc6e34')]:
            vv=[]
            for seed in range(200,204):
                z=[a for a in sr if a['mechanism']==mech and float(a['N2'])==N2 and int(a['seed'])==seed]
                t=[float(a['t']) for a in z]; vv.append([float(a['variance'])/.10125 for a in z])
            ax[j].plot(t,np.mean(vv,axis=0),label=label,color=color)
            ax[j].fill_between(t,np.min(vv,axis=0),np.max(vv,axis=0),alpha=.15,color=color)
        ax[j].set(title=f'{mech.capitalize()} buoyancy background',xlabel='Dimensionless time',ylabel='Solute variance / initial variance')
        ax[j].legend(fontsize=9)
    fig.suptitle('BOUSSINESQ SENSITIVITY · prescribed mean gradients, four seeds',fontsize=13)
    save(fig,'stratification')

    cr=read('freezing_model_curves.csv'); dr=read('freezing_simulated_droplets.csv')
    armspec=json.loads((ROOT/'data/freezing_arms.json').read_text())
    conditional=[]
    for row in cr:
        xx=float(row['cooling_below_reference_K'])
        if xx>=12: continue
        HH=hazard_model(np.array([xx,xx+.05]),**armspec[row['arm']])[1]
        conditional.append([row['arm'],xx,xx+.05,float(-np.expm1(-(HH[1]-HH[0])))])
    write_csv(ROOT/'data/freezing_conditional_probabilities.csv',['arm','cooling_start_K','cooling_end_K','P_freeze_in_interval_given_liquid_no_frailty'],conditional)
    fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    colors=['#718096','#167d9a','#cc6e34','#7449a3']
    arms=['distilled_matched_null','bio_added_assumed','mineral_added_assumed','bio_and_mineral_assumed']
    for arm,color in zip(arms,colors):
        q=[a for a in cr if a['arm']==arm]
        x=np.array([float(a['cooling_below_reference_K']) for a in q]); S=[float(a['survival']) for a in q]
        ax[0].plot(x,S,label=arm.replace('_assumed','').replace('_',' '),color=color)
        ax[1].semilogy(x,[float(a['hazard_per_min']) for a in q],color=color)
        ev=np.array([float(a['event_or_censor_cooling_K']) for a in dr if a['arm']==arm])
        ax[2].hist(ev,bins=np.linspace(0,12,61),histtype='step',density=True,color=color,label=arm)
    ax[0].legend(fontsize=7)
    for a in ax: a.set_xlabel('Cooling below unspecified reference (K)')
    ax[0].set(ylabel='Probability still liquid',title='Exact survival (no frailty)')
    ax[1].set(ylabel='Conditional event rate (min⁻¹)',ylim=(1e-5,1e4),title='Assumed component hazards')
    ax[2].set(ylabel='Simulated event density (K⁻¹)',title='12 batches × 96 drops, with frailty')
    fig.suptitle('ILLUSTRATIVE FREEZING MODEL · no absolute freezing-temperature predictions',fontsize=13)
    save(fig,'freezing')

    # Empirical survival & batch-bootstrap pointwise bands, including independent
    # null arms. Repeated droplets are nested within batches, never pooled as
    # independent experimental replicates for confidence intervals.
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    survival_rows=[]
    rng=np.random.default_rng(771)
    for j,selected in enumerate([list(json.loads((ROOT/'data/freezing_arms.json').read_text()))[:5],
                               ['bio_added_assumed','bio_slow_cooling','bio_fast_cooling']]):
        for arm in selected:
            x=np.linspace(0,12,121)
            bs=np.array([[np.mean(np.array([float(a['event_or_censor_cooling_K']) for a in dr if a['arm']==arm and int(a['batch'])==b])>v) for v in x] for b in range(12)])
            boot=bs[rng.integers(0,12,(2000,12))].mean(axis=1)
            lo,hi=np.quantile(boot,[.025,.975],axis=0); mean=bs.mean(axis=0)
            ax[j].plot(x,mean,label=arm.replace('_matched_null','').replace('_',' '))
            ax[j].fill_between(x,lo,hi,alpha=.12)
            survival_rows.extend([[arm,xx,mm,ll,hh] for xx,mm,ll,hh in zip(x,mean,lo,hi)])
        ax[j].set(xlabel='Cooling below unspecified reference (K)',ylabel='Simulated fraction liquid',title=['Treatment-label null','Cooling-rate sensitivity'][j])
        ax[j].legend(fontsize=8)
    ax[1].set_xlim(2,4)
    fig.suptitle('SIMULATED survival · pointwise 95% batch-bootstrap bands',fontsize=13)
    save(fig,'freezing_controls')
    write_csv(ROOT/'data/freezing_simulated_survival.csv',['arm','cooling_K','survival','bootstrap95_low','bootstrap95_high'],survival_rows)

    # Spatial solute/particle overlap on exactly the same 16x16 bins.
    overlaps=[]
    coarse=m['concentration'].reshape(-1,16,4,16,4).mean(axis=(2,4))
    # Histograms occupy [j/16,(j+1)/16), whereas spectral point samples are at
    # grid nodes. Block averaging is an approximate bin integral; identify it.
    for mode in ['inert','brownian','motile_assumed']:
        pp=np.load(ROOT/f'data/particles_{mode}_100.npz')['counts']/8192
        sol=coarse/coarse.sum(axis=(1,2))[:,None,None]
        val=1-.5*np.sum(abs(pp-sol),axis=(1,2))
        overlaps.extend([[mode,float(t),float(v)] for t,v in zip(m['time'],val)])
    write_csv(ROOT/'data/solute_particle_overlap_seed100.csv',['mode','t','overlap_approximate_common_bins'],overlaps)

    # Saved-data numerical quality checks; fail rather than silently publish NaN.
    numeric=[]
    for path in (ROOT/'data').glob('*.npz'):
        with np.load(path) as z:
            assert all(np.all(np.isfinite(z[k])) for k in z.files), path
    conv=read('stratification_convergence.csv')
    maxdiv=max(float(a['max_divergence']) for a in sr)
    minc=min(float(a['minimum']) for a in sr); maxc=max(float(a['maximum']) for a in sr)
    checks={'all_saved_arrays_finite':True,'stratification_max_divergence':maxdiv,
            'stratification_min_concentration':minc,'stratification_max_concentration':maxc,
            'stratification_max_mass_error':max(abs(float(a['mass'])-.5) for a in sr),
            'stratification_max_primary_low_mode_error':max(float(a['low_mode_L2_all_fields_vs_reference']) for a in conv if a['is_primary']=='True'),
            'stratification_max_finest_step_low_mode_error':max(float(a['low_mode_L2_all_fields_vs_reference']) for a in conv if a['is_finest_step']=='True')}
    assert maxdiv<1e-10 and minc>-1e-6 and maxc<1+1e-6
    assert checks['stratification_max_primary_low_mode_error']<.005
    assert checks['stratification_max_finest_step_low_mode_error']<.0005
    write_json(ROOT/'verification.json',checks)
    print(json.dumps(checks,indent=2))


if __name__=='__main__': main()
