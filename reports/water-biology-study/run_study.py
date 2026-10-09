"""Run from any directory: python run_study.py. Writes only within --out."""
from pathlib import Path
import argparse, csv, hashlib, json, platform, time
import numpy as np
from core import (scalar_run, particle_run, column_run, Boussinesq,
                  hazard_model, bootstrap_mean)

ROOT=Path(__file__).resolve().parent


def write_csv(path,header,rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(header); w.writerows(rows)


def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def run(out):
    started=time.perf_counter()
    data=out/'data'; data.mkdir(exist_ok=True,parents=True)
    configs={'status':'SIMULATIONS AND ANALYTIC CALCULATIONS; no new empirical samples',
      'solute':{'n':64,'dt':.01,'duration':8.,'D':.001,'seeds':list(range(100,108)),
                'domain':'unit 2D torus','initial':'.5+.45*cos(2*pi*x)',
                'flow':'alternating sine shears, switch every 0.5, 16 seeded phases',
                'mapping':{'L_m':.001,'U_m_s':.001,'time_unit_s':1.,'D_m2_s':1e-9}},
      'particles':{'count':8192,'seeds':list(range(100,108)),'radius_m':.5e-6,'temperature_K':293.15,
                   'mu_Pa_s':.001,'excess_density_kg_m3':50,'histogram_bins_per_axis':16,
                   'motility_assumed_m_s':20e-6,'rotational_diffusivity_assumed_s_inverse':1},
      'stratification':{'n':32,'dt':.01,'duration':4.,'nu':.01,'D':.001,
                         'density_kappa':.001,'thermal_kappa':.14,'N2':[-4,0,4],
                         'seeds':list(range(200,204)),
                         'unstable_density_override':{'n':64,'dt':.005},
                         'units':'nondimensional sensitivity study; not same water mapping as Q1'},
      'freezing':{'batches':12,'droplets_per_batch':96,'seeds':list(range(500,512)),
                  'volume_reference_uL':1.,'cooling_reference_K_min':1.,
                  'axis':'K below unspecified reference temperature; no absolute Celsius prediction',
                  'lognormal_batch_log_sd':.25,'lognormal_droplet_log_sd':.35,
                  'background_slope_assumed':1.5,'background_anchor_assumed':8.,
                  'mineral_slope_assumed':2.,'mineral_anchor_assumed':6.,
                  'bio_slope_literature_informed':8.7,'bio_anchor_assumed':3.}}
    write_json(out/'protocol.json',configs)
    scalar_rows=[]; scalar_summ=[]; primary={}
    for flow in ['still','stirred']:
        for seed in range(100,108):
            r,m=scalar_run(seed=seed,flow=flow)
            scalar_rows.extend([[flow,seed,*a] for a in r])
            idx=np.flatnonzero(r[:,2]/r[0,2]<=.05)
            t95=float(r[idx[0],0]) if len(idx) else None
            scalar_summ.append({'flow':flow,'seed':seed,'mixing':float(r[-1,3]),
                                'overlap':float(r[-1,4]),'t95':t95,
                                'mass_error':float(np.max(abs(r[:,1]-.5))),
                                'min':float(r[:,5].min()),'max':float(r[:,6].max())})
            if seed==100:
                np.savez_compressed(data/f'solute_{flow}_maps.npz',time=r[:,0],concentration=m)
                primary[flow]=m
    write_csv(data/'solute_timeseries.csv',['flow','seed','t','mass','variance','mixing_index','water_overlap','minimum','maximum'],scalar_rows)
    write_json(data/'solute_summary.json',scalar_summ)
    print('Solute ensembles complete',flush=True)

    conv=[]; fields={}
    for n,dt in [(32,.01),(64,.01),(128,.01),(128,.02),(128,.005)]:
        r,m=scalar_run(n=n,dt=dt,save_every=8.)
        fields[n,dt]=m[-1]
    for n,dt in fields:
        target=fields[128,.005]
        sampled=target[::128//n,::128//n]
        conv.append([n,dt,float(np.sqrt(np.mean((fields[n,dt]-sampled)**2))),float(fields[n,dt].var())])
    grid32=np.sqrt(np.mean((fields[32,.01]-fields[64,.01][::2,::2])**2))
    grid64=np.sqrt(np.mean((fields[64,.01]-fields[128,.01][::2,::2])**2))
    step20=np.sqrt(np.mean((fields[128,.02]-fields[128,.01])**2))
    step10=np.sqrt(np.mean((fields[128,.01]-fields[128,.005])**2))
    write_csv(data/'scalar_convergence.csv',['n','dt','rms_vs_128_dt005','final_variance'],conv)

    kb=1.380649e-23; T=293.15; mu=.001; radius=.5e-6
    Db=kb*T/(6*np.pi*mu*radius)
    particle_rows=[]; particle_summ=[]; count_checks=[]
    for mode in ['inert','brownian','motile_assumed']:
        for seed in range(100,108):
            r,h,x,e=particle_run(seed=seed,brownian=Db/1e-6 if mode!='inert' else 0.,
                                 swimming=.02 if mode=='motile_assumed' else 0.)
            particle_rows.extend([[mode,seed,*a] for a in r])
            particle_summ.append([mode,seed,*r[-1,1:]])
            np.savez_compressed(data/f'particles_{mode}_{seed}.npz',time=r[:,0],counts=h,
                                initial=x[0],final=x[-1],shear_exposure=e)
    for count in [2048,8192,32768]:
        for seed in range(100,108):
            r,_,_,_=particle_run(seed=seed,count=count,dt=.02)
            count_checks.append([count,seed,r[-1,2],r[-1,3]])
    write_csv(data/'particle_timeseries.csv',['mode','seed','t','cv2','cv2_minus_uniform_noise','overlap_with_uniform','mean_shear_exposure','p95_shear_exposure','count'],particle_rows)
    write_csv(data/'particle_summary.csv',['mode','seed','cv2','cv2_minus_uniform_noise','overlap_with_uniform','mean_shear_exposure','p95_shear_exposure','count'],particle_summ)
    write_csv(data/'particle_count_check.csv',['count','seed','corrected_cv2','overlap_with_uniform'],count_checks)
    print('Particle ensembles complete',flush=True)

    column=[]
    for mode,r,k in [('single',radius,0.),('larger_inert_proxy',5e-6,0.),('adsorption_assumed',radius,1/3600)]:
        for seed in range(300,316):
            velocity,budget,z0=column_run(seed,radius=r,adsorption_rate=k)
            column.append([mode,seed,r,k,velocity,*budget])
    write_csv(data/'column_budgets.csv',['mode','seed','radius_m','adsorption_rate_s_inverse','settling_m_s','suspended','deposited','attached'],column)

    strat=[]; sf={}
    for mechanism,kappa in [('density',.001),('thermal',.14)]:
        for N2 in [0.,4.,-4.]:
            for seed in range(200,204):
                n,dt=(64,.005) if mechanism=='density' and N2<0 else (32,.01)
                model=Boussinesq(n,kappa=kappa,N2=N2)
                r,f=model.run(seed,dt=dt)
                strat.extend([[mechanism,N2,seed,*a] for a in r])
                if seed==200:
                    sf[mechanism,N2,n,dt]=f[-1]
                    np.savez_compressed(data/f'stratification_{mechanism}_{N2:g}.npz',time=r[:,0],fields=f)
    stratconv=[]
    for mechanism,kappa in [('density',.001),('thermal',.14)]:
        for N2 in [4.,-4.]:
            unstable=mechanism=='density' and N2<0
            runs=[(32,.01),(48,.01),(64,.0025),(96,.0025),(128,.0025),(128,.00125)] if unstable else [(48,.01),(64,.005),(64,.0025)]
            for n,dt in runs:
                r,f=Boussinesq(n,kappa=kappa,N2=N2).run(dt=dt)
                sf[mechanism,N2,n,dt]=f[-1]
            ref=sf[mechanism,N2,128,.00125] if unstable else sf[mechanism,N2,64,.0025]
            checks=[(32,.01),(48,.01),(64,.005),(64,.0025),(96,.0025),(128,.0025)] if unstable else [(32,.01),(48,.01),(64,.005)]
            for n,dt in checks:
                f=sf[mechanism,N2,n,dt]
                # Compare low Fourier modes on a common 32-grid for all n.
                m=np.fft.fftfreq(32)*32
                def low(a):
                    inds=(m.astype(int)%a.shape[-1])
                    return np.fft.fft2(a)[:,inds[:,None],inds[None,:]]/a.shape[-1]**2
                err=float(np.sqrt(np.sum(abs(low(f)-low(ref))**2)))
                is_primary=(n,dt)==((64,.005) if unstable else (32,.01))
                is_finest_step=(n,dt)==((128,.0025) if unstable else (64,.005))
                stratconv.append([mechanism,N2,n,dt,err,float(f[3].var()),is_primary,is_finest_step])
    write_csv(data/'stratification_timeseries.csv',['mechanism','N2','seed','t','mass','variance','kinetic','stable_potential','max_divergence','minimum','maximum'],strat)
    write_csv(data/'stratification_convergence.csv',['mechanism','N2','n','dt','low_mode_L2_all_fields_vs_reference','variance','is_primary','is_finest_step'],stratconv)
    print('Stratification controls complete',flush=True)

    # Q2: dilute passive spheres and controlled material-property sensitivities.
    stress=[]
    for per_ml in [0,1e4,1e6,1e8]:
        phi=per_ml*1e6*4*np.pi*radius**3/3
        mur=1+2.5*phi
        stress.append([per_ml,phi,mur,mur,1/mur])
    write_csv(data/'passive_rheology.csv',['particles_per_ml','volume_fraction','mu_relative','stress_relative_at_fixed_shear','shear_relative_at_fixed_stress'],stress)
    thermal=[]
    for label,mur,sigmar,heatcapr,Kevapr in [('matched',1,1,1,1),('viscosity_only',1.01,1,1,1),
        ('surface_tension_only',1,.9,1,1),('heat_capacity_only',1,1,1.01,1),('evaporation_K_only',1,1,1,1.1)]:
        thermal.append([label,mur,sigmar,heatcapr,Kevapr,mur/sigmar,
                        1/np.sqrt(sigmar),heatcapr,1/Kevapr])
    write_csv(data/'stress_property_ablation.csv',['label','mu_relative','sigma_relative','heat_capacity_relative','D2law_K_relative','capillary_number_relative','inertial_capillary_time_relative','lumped_thermal_time_relative','D2law_lifetime_relative'],thermal)

    # Q3: no treatment-specific empirical rate is fabricated. First five arms
    # intentionally share identical physical hazards; sample draws are independent.
    arms={name:{} for name in ['distilled_matched_null','filtered_matched_null','heated_matched_null','boiled_matched_null','untreated_matched_null']}
    arms.update({'bio_added_assumed':{'bio':1.},'mineral_added_assumed':{'mineral':1.},
       'bio_and_mineral_assumed':{'bio':1.,'mineral':1.},'bio_removed_ablation':{},
       'surface_x10_assumed':{'surface':10.},'bio_slow_cooling':{'bio':1.,'cooling':.5},
       'bio_fast_cooling':{'bio':1.,'cooling':2.},'bio_half_volume':{'bio':1.,'volume':.5},
       'bio_double_volume':{'bio':1.,'volume':2.}})
    x=np.linspace(0,12,2401); curves=[]; samples=[]; batches=[]
    for ai,(arm,kw) in enumerate(arms.items()):
        rate,H,S=hazard_model(x,**kw)
        for xx,rr,hh,ss in zip(x[::10],rate[::10],H[::10],S[::10]):
            curves.append([arm,xx,rr,hh,ss])
        for batch in range(12):
            rng=np.random.default_rng(500+batch+1000*ai)
            bm=rng.lognormal(-.25**2/2,.25)
            frailty=bm*rng.lognormal(-.35**2/2,.35,96)
            threshold=rng.exponential(size=96)/frailty
            events=np.interp(threshold,H,x)
            observed=threshold<=H[-1]
            for j in range(96):
                samples.append([arm,batch,j,bm,frailty[j],events[j],int(observed[j])])
            batches.append([arm,batch,float(np.mean(events<=5)),float(np.median(events)),int(observed.sum())])
    write_csv(data/'freezing_model_curves.csv',['arm','cooling_below_reference_K','hazard_per_min','cumulative_hazard','survival'],curves)
    write_csv(data/'freezing_simulated_droplets.csv',['arm','batch','droplet','batch_multiplier','total_rate_multiplier','event_or_censor_cooling_K','event_observed'],samples)
    write_csv(data/'freezing_simulated_batches.csv',['arm','batch','fraction_frozen_by_5K','median_event_cooling_K','events'],batches)
    write_json(data/'freezing_arms.json',arms)
    # Published slope uncertainty propagated as a sensitivity band, NOT new CI.
    shifts=[]
    for ratio in [.1,.5,1,2,10,100]:
        shifts.append([ratio,-np.log(ratio)/8.7,*sorted([-np.log(ratio)/8.3,-np.log(ratio)/9.1])])
    write_csv(data/'published_slope_cooling_shifts.csv',['cooling_rate_ratio','predicted_delta_T50_K','sensitivity_low_K','sensitivity_high_K'],shifts)
    print('Freezing probability controls complete',flush=True)

    summary={'solute':{},'particles':{},'stratification':{},'freezing_simulated':{},
      'numerics':{'scalar_grid32_to64_L2':float(grid32),'scalar_grid64_to128_L2':float(grid64),
       'scalar_dt02_to01_L2':float(step20),'scalar_dt01_to005_L2':float(step10),
       'temporal_observed_order':float(np.log2(step20/step10)),
       'max_solute_mass_error':max(r['mass_error'] for r in scalar_summ)},
      'physics':{'brownian_D_m2_s':float(Db),'solute_to_cell_diffusion_ratio':float(1e-9/Db),
       'cell_1d_rms_displacement_8s_um':float(np.sqrt(2*Db*8)*1e6),
       'solute_1d_rms_displacement_8s_um':float(np.sqrt(2e-9*8)*1e6),
       'single_cell_settling_m_s':float(column[0][4]),
       'passive_viscosity_increase_at_1e6_per_ml_ppm':float((stress[2][2]-1)*1e6)}}
    for flow in ['still','stirred']:
        ss=[a for a in scalar_summ if a['flow']==flow]
        summary['solute'][flow]={'mixing_mean_bootstrap95':bootstrap_mean([a['mixing'] for a in ss]),
             't95': [a['t95'] for a in ss]}
    for mode in ['inert','brownian','motile_assumed']:
        rr=[a for a in particle_summ if a[0]==mode]
        summary['particles'][mode]={'corrected_cv2_mean_bootstrap95':bootstrap_mean([a[3] for a in rr]),
                                     'overlap_mean_bootstrap95':bootstrap_mean([a[4] for a in rr])}
    for mechanism in ['density','thermal']:
        for N2 in [0.,4.,-4.]:
            rr=[a for a in strat if a[0]==mechanism and a[1]==N2 and a[3]==4.]
            summary['stratification'][f'{mechanism}_{N2:g}']={'final_variance_mean_bootstrap95':bootstrap_mean([a[5] for a in rr])}
    for arm in arms:
        rr=[a for a in batches if a[0]==arm]
        summary['freezing_simulated'][arm]={'fraction_by_5K_mean_batch_bootstrap95':bootstrap_mean([a[2] for a in rr]),
                                           'batch_median_cooling_mean_bootstrap95':bootstrap_mean([a[3] for a in rr])}
    write_json(out/'summary.json',summary)
    import scipy, matplotlib, pytest
    write_json(out/'environment.json',{'python':platform.python_version(),'numpy':np.__version__,
        'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'pytest':pytest.__version__,
        'platform':platform.platform(),'elapsed_seconds':time.perf_counter()-started,
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}})
    print(json.dumps(summary['numerics'],indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--out',type=Path,default=ROOT)
    run(parser.parse_args().out)
