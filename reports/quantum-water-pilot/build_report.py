"""Regenerate figures, readable report and summary tables from recorded data."""
import argparse,base64,csv,hashlib,html,json,platform,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import c,h,Avogadro
from vendor.ns_solver import Solver,difference

COLORS={'H2O':'#176d9c','D2O':'#cf6b2f'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':'#263b4d','text.color':'#263b4d','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','grid.alpha':.2})

def read(p):return json.loads(p.read_text())
def dump(p,d):p.write_text(json.dumps(d,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()),encoding='utf-8')
def csvfile(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def savefig(fig,root,name):
    fig.savefig(root/f'{name}.png',dpi=170,bbox_inches='tight');fig.savefig(root/f'{name}.svg',bbox_inches='tight');plt.close(fig)
def table(headers,rows):
    return '<div class="table"><table><thead><tr>'+''.join('<th>'+html.escape(str(s))+'</th>' for s in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(s))+'</td>' for s in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def img(root,name,caption):
    encoded=base64.b64encode((root/f'{name}.png').read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{encoded}" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>'

def vibration(root):
    cfg=read(root/'protocol.json');records=read(root/'vibration-final/runs.json');mode=cfg['mode_order'];rows=[];pred={};obs={}
    for iso in COLORS:
        r=next(x for x in records if x['method']=='GFN2-xTB' and x['isotope']==iso and x['hessian_step_bohr']==.00125)
        pred[iso]=np.array(r['frequency_cm1']);obs[iso]=np.array(cfg[iso+'_observed_cm-1'])
        for j,m in enumerate(mode):rows.append(dict(isotope=iso,mode=m,predicted_cm1=pred[iso][j],observed_IR_cm1=obs[iso][j],error_cm1=pred[iso][j]-obs[iso][j],relative_error_percent=100*(pred[iso][j]/obs[iso][j]-1)))
    shifts=[]
    oxygen=cfg['mass_u']['O16'];hydrogen=cfg['mass_u']['H1'];deuterium=cfg['mass_u']['D2']
    ratio=np.sqrt((oxygen*hydrogen/(oxygen+hydrogen))/(oxygen*deuterium/(oxygen+deuterium)))
    for j,m in enumerate(mode):
        shifts.append(dict(mode=m,predicted_H_minus_D_cm1=pred['H2O'][j]-pred['D2O'][j],observed_H_minus_D_cm1=obs['H2O'][j]-obs['D2O'][j],predicted_redshift_percent=100*(1-pred['D2O'][j]/pred['H2O'][j]),observed_redshift_percent=100*(1-obs['D2O'][j]/obs['H2O'][j]),classical_local_oscillator_redshift_percent=100*(1-ratio)))
    errors=np.array([r['error_cm1'] for r in rows]);relative=np.array([r['relative_error_percent'] for r in rows])
    stats=dict(RMSE_cm1=float(np.sqrt(np.mean(errors**2))),MAE_cm1=float(np.mean(abs(errors))),MAPE_percent=float(np.mean(abs(relative))),max_relative_error_percent=float(np.max(abs(relative))),primary_gate=read(root/'physics-checks.json')['vibration_gate'],quantum_harmonic_ZPE_kJmol={i:float(.5*np.sum(pred[i])*h*c*100*Avogadro/1000) for i in pred})
    secondary=[r for r in records if r['method']=='GFN1-xTB']
    stats['GFN1_max_error_percent']=max(float(np.max(abs(np.array(r['frequency_cm1'])/obs[r['isotope']]-1)))*100 for r in secondary)
    href=np.array([1649.,3832.,3943.])
    stats['H2O_harmonic_reference_cm1']=href
    stats['H2O_harmonic_reference_error_percent']=100*(pred['H2O']/href-1)
    stats['H2O_harmonic_reference_max_error_percent']=float(np.max(abs(pred['H2O']/href-1))*100)
    csvfile(root/'vibration_comparison.csv',rows);csvfile(root/'isotope_shifts.csv',shifts);dump(root/'vibration-summary.json',stats)
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    ax=axes[0];x=np.array([1000,4000]);ax.fill_between(x,.95*x,1.05*x,color='#edf1f3',label='±5% pilot screen');ax.plot(x,x,color='#75828a',lw=1)
    for iso in COLORS:ax.scatter(obs[iso],pred[iso],s=55,label=iso,color=COLORS[iso])
    ax.set(xlabel='Published IR position (cm⁻¹)',ylabel='Predicted harmonic position (cm⁻¹)',title='Six gas-phase band positions');ax.legend(fontsize=8)
    ax=axes[1];pos=np.arange(3);w=.26
    for offset,key,label,color in [(-w,'observed_redshift_percent','Observed','#263b4d'),(0,'predicted_redshift_percent','GFN2 + isotope masses','#176d9c'),(w,'classical_local_oscillator_redshift_percent','Local oscillator mass law','#bcb5a1')]:ax.bar(pos+offset,[r[key] for r in shifts],w,label=label,color=color)
    ax.set(xticks=pos,xticklabels=['Bend','Sym. stretch','Anti. stretch'],ylabel='H → D frequency decrease (%)',title='An isotope shift is not uniquely quantum',ylim=(0,34));ax.legend(fontsize=7)
    ax=axes[2]
    for j,r in enumerate(rows):ax.barh(j,r['error_cm1'],color=COLORS[r['isotope']])
    ax.axvline(0,color='black',lw=.6);ax.set(yticks=np.arange(6),yticklabels=[r['isotope']+' '+r['mode'].replace('symmetric','sym.').replace('antisym.','antisym.') for r in rows],xlabel='Prediction minus measurement (cm⁻¹)',title='Model error dominates numerical error')
    savefig(fig,root,'vibration_comparison')
    fig,axes=plt.subplots(2,1,figsize=(10,5),sharex=True,layout='constrained')
    for ax,iso in zip(axes,COLORS):
        ax.vlines(obs[iso],0,.9,color='#263b4d',lw=2,label='Published gas IR band positions');ax.vlines(pred[iso],0,-.9,color=COLORS[iso],lw=2,label='Calculated harmonic modes')
        ax.set(yticks=[],ylim=(-1.1,1.1),title=iso);ax.axhline(0,color='#a0a0a0',lw=.7);ax.legend(fontsize=8,loc='center right')
    axes[-1].set_xlabel('Wavenumber (cm⁻¹)');savefig(fig,root,'IR_band_positions')
    return rows,shifts,stats

def evidence_diagram(root):
    from matplotlib.patches import FancyBboxPatch
    fig,ax=plt.subplots(figsize=(12,4.6),layout='constrained');ax.set(xlim=(0,1),ylim=(-.10,1));ax.axis('off')
    def box(x,y,text,color='#eaf3f7',width=.27):
        ax.add_patch(FancyBboxPatch((x-width/2,y-.10),width,.20,boxstyle='round,pad=0.012',facecolor=color,edgecolor='#a3b7c2'))
        ax.text(x,y,text,ha='center',va='center',fontsize=11)
    def arrow(a,b,dashed=False):
        ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',lw=1.8,color='#b3582b' if dashed else '#176d9c',linestyle='--' if dashed else '-'))
    box(.2,.83,'Quantum electronic model\nGFN2-xTB')
    box(.7,.83,'Gas infrared band positions\nMeasured quantity predicted',width=.36)
    arrow((.35,.83),(.51,.83));ax.text(.43,.91,'calculated',ha='center',fontsize=9)
    box(.16,.46,'Liquid molecular model\nq-TIP4P/F + PIMD')
    box(.50,.46,'Viscosity from molecules\nNot validated here','#fff0e7')
    box(.84,.46,'Navier–Stokes flow\nVelocity and dissipation')
    arrow((.31,.46),(.35,.46),True);arrow((.65,.46),(.69,.46),True)
    box(.50,.12,'Experimental property\ncorrelations: density, viscosity')
    arrow((.65,.12),(.84,.34));ax.text(.78,.14,'actual inputs',ha='center',fontsize=9,color='#176d9c')
    ax.text(.02,-.075,'Solid arrows: completed connection. Dashed arrows: missing predictive link.',fontsize=9)
    savefig(fig,root,'evidence_chain')

def flow(root):
    p=root/'fluid-data';props=read(p/'properties.json');runs=read(p/'channel_results.json');conv=read(p/'channel_convergence.json');vr=read(p/'vortex_results.json');bench=read(p/'historical_benchmark.json');extended=read(p/'extended_convergence.json')
    real={r['material']:r for r in runs if r['T_K']==298.15 and r['G_Pam']==10 and r['material'] in COLORS}
    ratio=props['298.15']['H2O']['mu']/props['298.15']['D2O']['mu']
    fields=[];pert=[]
    for T in [283.15,298.15,313.15]:
        for U in [.005,.01,.02]:
            n=32 if U==.02 else 16
            if n==32:
                h=np.load(p/f'extended-T{T}-H2O.npz')['n32'];d=np.load(p/f'extended-T{T}-D2O.npz')['n32']
            else:
                h=np.load(p/f'vortex-T{T}-H2O-U{U}-n16.npz')['final_hat'];d=np.load(p/f'vortex-T{T}-D2O-U{U}-n16.npz')['final_hat']
            weight=np.full((1,1,h.shape[-1]),2.);weight[...,0]=weight[...,-1]=1.
            norm=lambda a:np.sqrt(np.sum(abs(a)**2*weight)/n**6)
            fields.append(dict(T_K=T,U_scale_ms=U,grid_n=n,rms_velocity_difference_ms=U*norm(d-h),relative_velocity_field_difference=norm(d-h)/norm(h)))
    for iso in COLORS:
        base=np.load(p/f'vortex-T298.15-{iso}-U0.01-n24-dt0.01-s0.npz');weight=np.full((1,1,13),2.);weight[...,0]=weight[...,-1]=1.
        norm=lambda a:np.sqrt(np.sum(abs(a)**2*weight)/24**6)
        for seed in [1,2,3]:
            other=np.load(p/f'vortex-T298.15-{iso}-U0.01-n24-dt0.01-s{seed}.npz')
            gain=norm(other['final_hat']-base['final_hat'])/norm(other['initial_hat']-base['initial_hat'])
            decay=norm(base['final_hat'])/norm(base['initial_hat'])
            pert.append(dict(isotope=iso,seed=seed,absolute_perturbation_amplitude_gain=gain,perturbation_energy_gain=gain**2,base_amplitude_decay=decay,relative_to_base_amplitude_gain=gain/decay))
    stats=dict(steady_channel_D_over_H=ratio,steady_channel_velocity_reduction_percent=100*(1-ratio),tau_D_over_H=real['D2O']['tau_s']/real['H2O']['tau_s'],channel_worst_relative_error=max(r['relative_L2_error'] for r in runs),channel_worst_energy_residual=max(r['energy_relative_residual'] for r in runs),vortex_worst_divergence=max(r['divergence_rms'] for r in vr),vortex_worst_energy_residual=max(abs(r['energy_residual']) for r in vr),strongest_forcing_worst_fine_field_difference=max(r['n24_n32'] for r in extended),historical_max_viscosity_ratio_error_percent=max(abs(r['relative_error_percent']) for r in bench),vortex_velocity_differences=fields,perturbation_controls=pert)
    dump(root/'fluid-summary.json',stats);csvfile(root/'channel_results.csv',runs);csvfile(root/'vortex_velocity_differences.csv',fields);csvfile(root/'historical_viscosity_comparison.csv',bench)
    fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    for iso in COLORS:
        data=np.load(p/f'channel-T298.15-{iso}-G10.0.npz')
        axes[0,0].plot(data['times_s'],data['trace'][:,0]*1000,color=COLORS[iso],label=iso)
        axes[0,1].plot(data['profiles_ms'][-1]*1000,data['y_m']*1000,color=COLORS[iso],label=iso)
    axes[0,0].set(xlabel='Time (s)',ylabel='Mean velocity (mm/s)',title='Matched channel startup, 25 °C');axes[0,0].legend()
    axes[0,1].set(xlabel='Velocity (mm/s)',ylabel='Across 1 mm gap (mm)',title='Velocity profile after 1 s')
    temp=[10,25,40];ratios=[props[str(T+273.15)]['H2O']['mu']/props[str(T+273.15)]['D2O']['mu'] for T in temp]
    axes[0,2].plot(temp,ratios,'o-',color=COLORS['D2O']);axes[0,2].set(xlabel='Temperature (°C)',ylabel='Steady velocity D₂O / H₂O',title='Difference persists across temperature',ylim=(.74,.88))
    labels=['H2O','density_only','viscosity_only','D2O']
    values=[next(r['rms_velocity_ms'] for r in vr if r['T_K']==298.15 and r['U_ms']==.01 and r['n']==16 and r['material']==i)*1000 for i in labels]
    axes[1,0].bar(np.arange(4),values,color=['#176d9c','#6baecb','#dda979','#cf6b2f']);axes[1,0].set(xticks=np.arange(4),xticklabels=['H₂O','D density','D viscosity','D₂O'],ylabel='RMS velocity (mm/s)',title='3D vortex: separate material properties')
    axes[1,1].plot([r['T_C'] for r in bench],[r['observed_viscosity_ratio'] for r in bench],'o',color='#263b4d',label='1949 capillary measurements')
    axes[1,1].plot([r['T_C'] for r in bench],[r['predicted_viscosity_ratio'] for r in bench],'x--',color='#cf6b2f',label='IAPWS correlation');axes[1,1].set(xlabel='Temperature (°C)',ylabel='Dynamic viscosity D₂O / H₂O',title='Constitutive consistency, not holdout validation');axes[1,1].legend(fontsize=8)
    for iso in COLORS:
        xs=[r for r in conv if r['material']==iso and r['dt_s']==.001]
        axes[1,2].loglog([r['n']-1 for r in xs],[r['relative_L2_error'] for r in xs],'o-',label=iso,color=COLORS[iso])
    axes[1,2].set(xlabel='Across-gap intervals',ylabel='Relative error against exact startup',title='Channel grid convergence at 0.1 s');axes[1,2].legend()
    savefig(fig,root,'fluid_connection')
    return stats,real,bench

def molecular(root):
    data=read(root/'molecular-data/summary.json');groups=data['groups'];runs=data['runs']
    assert data['complete_main_seeds']==[1,2,3], 'Do not publish an incomplete ensemble'
    main=[r for r in runs if r['n']==64 and r['beads'] in [1,32] and r['T_K']==300 and r['dt_ps']==.0005 and r['mobility_NVE_ps']==5]
    columns=['tag','T_K','n','beads','seed','bond_mean_A','intramolecular_RDF_peak_A','hydrogen_bond_degree','tetrahedral_order','coordination_035nm','nearest_OO_A']
    csvfile(root/'molecular_results.csv',[{k:r[k] for k in columns} for r in runs])
    csvfile(root/'mobility_diagnostics.csv',[{k:r[k] for k in ['tag','short_D_1_2ps_1e9_m2s','short_D_2_4ps_1e9_m2s','diffusion_window_relative_disagreement','diffusion_window_pass','NVE_energy_drift_NkBT','NVE_energy_drift_pass','neighbor_replaced_fraction_1ps']} for r in runs if 'short_D_1_2ps_1e9_m2s' in r])
    stats=dict(completed_runs=len(runs),completed_main_runs=len(main),main_window_passes=sum(r['diffusion_window_pass'] for r in main),main_energy_passes=sum(r['NVE_energy_drift_pass'] for r in main),structural_controls_passes=sum(r['structural_screen_pass'] for r in data['controls']),structural_controls_total=len(data['controls']),total_simulated_ps=sum(r['equilibration_ps']+r['production_ps']+r['mobility_NVE_ps'] for r in runs),bulk_diffusion_validated=False,viscosity_computed=False)
    dump(root/'molecular-summary.json',stats)
    fig,ax=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    keys=['bond_mean_A','hydrogen_bond_degree','tetrahedral_order'];labels=['Mean covalent bond (Å)','Hydrogen-bond degree','Tetrahedral order']
    for col,(key,label) in enumerate(zip(keys,labels)):
        for j,iso in enumerate(COLORS):
            vals=[groups[f'{iso}-P{b}'][key] for b in [1,32]]
            ax[0,col].errorbar(np.arange(2)+j*.12-.06,[v['mean'] for v in vals],yerr=[(v['CI95'][1]-v['CI95'][0])/2 for v in vals],fmt='o-',capsize=4,color=COLORS[iso],label=iso)
        ax[0,col].set(xticks=[0,1],xticklabels=['Classical','32-bead quantum'],ylabel=label,title=label+' at 300 K');ax[0,col].legend(fontsize=8)
    for iso in COLORS:
        for beads,style in [(1,'--'),(32,'-')]:
            cases=[r for r in main if r['isotope']==iso and r['beads']==beads];raw=[np.load(root/'molecular-data'/r['tag']/'analysis.npz') for r in cases]
            ax[1,0].plot(raw[0]['rdf_r_nm']*10,np.mean([r['rdf_gOO'] for r in raw],axis=0),style,color=COLORS[iso],label=f'{iso}, P={beads}')
    ax[1,0].set(xlabel='O–O distance (Å)',ylabel='gOO(r)',title='Local organization; 64-molecule box',xlim=(2.2,5.));ax[1,0].legend(fontsize=8)
    for j,iso in enumerate(COLORS):
        vals=[groups[f'{iso}-P{b}']['intramolecular_RDF_peak_A'] for b in [1,32]]
        ax[1,1].errorbar(np.arange(2)+j*.12-.06,[v['mean'] for v in vals],yerr=[(v['CI95'][1]-v['CI95'][0])/2 for v in vals],fmt='o-',color=COLORS[iso],label=f'{iso} calculated')
        ref=.990 if iso=='H2O' else .985;ax[1,1].axhspan(ref-.005,ref+.005,color=COLORS[iso],alpha=.10);ax[1,1].axhline(ref,color=COLORS[iso],ls=':',lw=1)
    ax[1,1].set(xticks=[0,1],xticklabels=['Classical','32-bead quantum'],ylabel='Intramolecular radial peak (Å)',title='Approximate comparison to neutron distances');ax[1,1].legend(fontsize=8)
    for iso in COLORS:
        for beads,style in [(1,'--'),(32,'-')]:
            rows=[r for r in runs if r['isotope']==iso and r['beads']==beads and r['seed']==1 and r['n']==64 and r['dt_ps']==.0005]
            rows=sorted(rows,key=lambda r:r['T_K'])
            ax[1,2].plot([r['T_K'] for r in rows],[r['hydrogen_bond_degree'] for r in rows],style+'o',color=COLORS[iso],label=f'{iso}, P={beads}')
    ax[1,2].set(xlabel='Temperature (K)',ylabel='Hydrogen-bond degree',title='Temperature sensitivity: one run per condition');ax[1,2].legend(fontsize=8)
    savefig(fig,root,'molecular_structure')
    fig,ax=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    for r in main:
        style='o' if r['beads']==32 else 'x'
        ax[0].scatter(r['short_D_1_2ps_1e9_m2s'],r['short_D_2_4ps_1e9_m2s'],marker=style,color=COLORS[r['isotope']])
    lim=max(max(r['short_D_1_2ps_1e9_m2s'],r['short_D_2_4ps_1e9_m2s']) for r in main)*1.15
    ax[0].plot([0,lim],[0,lim],color='gray',lw=1);ax[0].set(xlabel='MSD fit at 1–2 ps (10⁻⁹ m²/s)',ylabel='MSD fit at 2–4 ps (10⁻⁹ m²/s)',title='Diffusion is not yet a converged observable',xlim=(0,lim),ylim=(0,lim));ax[0].text(.03,.95,'Blue H₂O; orange D₂O\n× classical; ● quantum',transform=ax[0].transAxes,va='top',fontsize=8)
    mobility=[r for r in runs if 'NVE_energy_drift_NkBT' in r]
    mlabels=[('H' if r['isotope']=='H2O' else 'D')+'\n'+('C' if r['beads']==1 else 'Q')+str(r['seed'])+('\n½ step' if r['dt_ps']==.00025 else '') for r in mobility]
    ax[1].bar(np.arange(len(mobility)),[abs(r['NVE_energy_drift_NkBT']) for r in mobility],color=[COLORS[r['isotope']] for r in mobility]);ax[1].axhline(.005,color='black',ls='--',label='Acceptance screen');ax[1].set(xticks=np.arange(len(mobility)),xticklabels=mlabels,xlabel='C: classical; Q: quantum; number: seed',ylabel='Absolute energy drift / NkBT',title='Energy check includes failures');ax[1].tick_params(axis='x',labelsize=7);ax[1].legend(fontsize=8)
    controls=data['controls'];names=[]
    for r in controls:
        cfg=next(s for s in runs if s['tag']==r['tag']);names.append(f"N={cfg['n']}, P={cfg['beads']}\nΔt={cfg['dt_ps']*1000:g} fs")
    ax[2].bar(np.arange(len(controls)),[r['hbond_difference'] for r in controls],color='#798c97');ax[2].axhline(.15,color='black',ls='--');ax[2].axhline(-.15,color='black',ls='--');ax[2].set(xticks=np.arange(len(controls)),xticklabels=names,ylabel='Change from main H₂O run, seed 1',title='Bead, step and size sensitivity');ax[2].tick_params(axis='x',labelsize=8)
    savefig(fig,root,'molecular_stress')
    return data,stats

SOURCES=[
('nist_h','H₂O gas IR, NIST WebBook / Shimanouchi 1972','https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=800'),
('nist_d','D₂O gas IR, NIST WebBook / Shimanouchi 1972','https://webbook.nist.gov/cgi/cbook.cgi?ID=C7789200&Mask=800'),
('nist_hdo','HDO gas IR, NIST WebBook / Shimanouchi 1972','https://webbook.nist.gov/cgi/cbook.cgi?ID=C14940637&Mask=800'),
('harmonic','NIST CCCBDB H₂O harmonic frequencies and equilibrium geometry','https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0'),
('xtb','xTB Hessian and isotope documentation','https://xtb-docs.readthedocs.io/en/latest/hessian.html'),
('gfn2','Bannwarth, Ehlert & Grimme (2019), GFN2-xTB quantum chemical method','https://doi.org/10.1021/acs.jctc.8b01176'),
('xtb_issue','xTB issue 1143: isotopic masses in numerical Hessian calculations','https://github.com/grimme-lab/xtb/issues/1143'),
('qtip','Habershon, Markland & Manolopoulos, JCP 131, 024501 (2009), q-TIP4P/F','https://arxiv.org/pdf/1011.1047'),
('neutron','Zeidler et al., JPCM 24, 284126 (2012), neutron isotope effects','https://purehost.bath.ac.uk/ws/portalfiles/portal/9301378/18O_JPCM_20111126.pdf'),
('rpmd','OpenMM RPMDIntegrator API','https://docs.openmm.org/latest/api-python/generated/openmm.openmm.RPMDIntegrator.html'),
('iapws','IAPWS releases: water and heavy-water material properties','https://iapws.org/documents/release'),
('visc','IAPWS 2020 heavy-water viscosity formulation, verification table 3','https://iapws.org/public/documents/OM7f5/D2Ovisc.pdf'),
('hardy','Hardy & Cottington, J. Research NBS 42, 573 (1949), table 2','https://nvlpubs.nist.gov/nistpubs/jres/42/jresv42n6p573_A1b.pdf'),
('visc_fit','Assael et al. (2021), heavy-water viscosity correlation and fitted datasets','https://srd.nist.gov/jpcrdreprint/5.0048711.pdf'),
('later_visc','Ragueneau, Caupin & Issenmann (2022), supplemental viscosity measurements','https://arxiv.org/pdf/2112.09024'),
('interface_flow','Yoon & Yoon (2018), silica-interface HPLC flow observations','https://www.nature.com/articles/s41598-018-24886-y'),
('mills','Mills, J. Physical Chemistry 77, 685 (1973), table III','https://fenix.ciencias.ulisboa.pt/downloadFile/2533412229372806/Mills1973.pdf'),
('repo','Existing project, inspected snapshot 4e4a970e20d3289e883dec5c34e1a1bd3968795e','https://github.com/NousVolition/My-Sources-Project-ALL/tree/4e4a970e20d3289e883dec5c34e1a1bd3968795e'),
]

def report(root,vrows,shifts,v,flowstats,real,bench,md,mstats):
    physics=read(root/'physics-checks.json');extra=read(root/'isotope-stress.json');groups=md['groups'];runs=md['runs'];sens=read(root/'fluid-data/sensitivity.json')
    stress=read(root/'fluid-stress/convergence.json');stressrows=read(root/'fluid-stress/results.json')
    assert len(stressrows)==12 and len(runs)==24, 'Missing a planned stress run'
    refinement=read(root/'fluid-refinement/convergence.json')
    assert len(read(root/'fluid-refinement/results.json'))==6, 'Missing a finer-grid follow-up'
    finest=[r for r in refinement if r['kind']=='grid' and r['fine_n']==96]
    primary=read(root/'vibration-final/runs.json');force=read(root/'molecular-force-validation.json')
    def ci(v,d=4):return f"{v['mean']:.{d}f} [{v['CI95'][0]:.{d}f}, {v['CI95'][1]:.{d}f}]" if v['CI95'] else f"{v['mean']:.{d}f}"
    def ref(id,label='source'):return f'<a href="#{id}">{label}</a>'
    content=[]
    content.append(f'''<header><p class="eyebrow">NOUSVOLITION • COMPUTATIONAL RESEARCH PILOT • OCTOBER 2026</p><h1>From quantum water<br>to observable fluid behavior</h1><p class="lede">The isolated-molecule prediction passes its pilot test. Measured isotope-dependent properties produce clear flow differences. The simulated molecular-to-fluid causal bridge remains unvalidated.</p></header>
    <div class="cards"><div><b>{v['max_relative_error_percent']:.2f}%</b><span>largest H₂O/D₂O IR-position error</span></div><div><b>{flowstats['steady_channel_velocity_reduction_percent']:.2f}%</b><span>slower D₂O steady channel flow at 25 °C</span></div><div><b>{mstats['completed_runs']}</b><span>molecular runs, including sensitivity controls</span></div></div>
    <h2>What the completed tests support</h2><p>Quantum electronic modeling predicts a measurable gas-phase isotope shift without fitting these results to the listed IR positions. On the liquid model, nuclear quantum sampling changes molecular geometry, but short trajectories and failed control screens prevent a claim of improved bulk dynamics. Continuum simulations quantify what measured density and viscosity differences do to flow. These are three separate calculations: viscosity was not derived from the quantum molecular simulations.</p>
    <p>These tests support separate parts of the central hypothesis; the complete causal chain remains unvalidated. This pilot does <strong>not</strong> establish that nuclear quantum effects are essential to the observed flow difference, or that one validated calculation now connects electronic structure to macroscopic transport.</p>''')
    ledger=[
        dict(test='Primary H₂O/D₂O IR positions',status='PASS',evidence=f"Maximum error {v['max_relative_error_percent']:.3f}% < 5%"),
        dict(test='Isotope shifts and numerical repeatability',status='PASS',evidence=f"{100*physics['vibration_gate']['max_isotope_shift_fraction_error']:.3f} percentage points; {physics['vibration_gate']['max_numerical_spread_cm1']:.3f} cm⁻¹ spread"),
        dict(test='HDO generalization without refitting',status='PASS',evidence=f"Maximum error {extra['max_absolute_error_percent']:.3f}% < 5%"),
        dict(test='Alternative GFN1 method',status='FAIL',evidence=f"Maximum error {v['GFN1_max_error_percent']:.3f}% > 5%"),
        dict(test='Like-for-like H₂O harmonic reference',status='FAIL',evidence=f"Maximum error {v['H2O_harmonic_reference_max_error_percent']:.3f}% > 5%; exposes error cancellation"),
        dict(test='Liquid structure controls',status='LIMITED',evidence=f"{mstats['structural_controls_passes']}/{mstats['structural_controls_total']} single-run screens pass; sampling bias unresolved"),
        dict(test='Bulk diffusion / viscosity from molecules',status='UNVALIDATED',evidence=f"Diffusion-window checks {mstats['main_window_passes']}/12; energy checks {mstats['main_energy_passes']}/12; no viscosity estimate"),
        dict(test='Channel analytic verification',status='PASS',evidence=f"Worst relative profile error {flowstats['channel_worst_relative_error']:.2g}"),
        dict(test='Original 3D range, refined grids',status='PASS',evidence=f"Worst final field change {100*flowstats['strongest_forcing_worst_fine_field_difference']:.3f}% < 1%"),
        dict(test='100 mm/s, original 32³ → 48³ comparison',status='FAIL',evidence=f"Field change up to {100*max(r['relative_field_change'] for r in stress if r['fine_n']==48 and r['U_ms']==.1):.3f}% > 1%"),
        dict(test='100 mm/s, follow-up 64³ → 96³ comparison',status='PASS' if all(r['passed'] for r in finest) else 'FAIL',evidence=f"Field change up to {100*max(r['relative_field_change'] for r in finest):.3f}%; same 1% screen")]
    dump(root/'stress_summary.json',ledger);content.append(table(['Test','Result','Evidence'],[[r['test'],r['status'],r['evidence']] for r in ledger]))
    evidence_diagram(root);content.append(img(root,'evidence_chain','The gas and liquid calculations use separate potentials. No computed molecular viscosity is passed to the continuum solver.'))
    content.append('<h2>1. Gas-phase vibration: prediction versus measurement</h2><p>GFN2-xTB is a semiempirical tight-binding quantum chemical method. It supplies an electronic potential and Cartesian Hessian at an optimized isolated-water geometry. Independent mass weighting gives the three internal harmonic modes for each isotope. The primary comparison uses the finest Hessian displacement, 0.00125 bohr, with no frequency scaling. The observations are published gas-phase infrared fundamental positions, not liquid absorption peaks. '+ref('gfn2','Method')+'; '+ref('nist_h','H₂O data')+'; '+ref('nist_d','D₂O data')+'.</p>')
    content.append(table(['Isotope','Mode','Calculated / cm⁻¹','Observed IR / cm⁻¹','Error / %'],[[r['isotope'],r['mode'],f"{r['predicted_cm1']:.2f}",f"{r['observed_IR_cm1']:.2f}",f"{r['relative_error_percent']:+.3f}"] for r in vrows]))
    content.append(img(root,'vibration_comparison','Unscaled predictions compared with experimental band positions. The simple local-oscillator mass law is only a diagnostic baseline; the same full Hessian with classical nuclei reproduces the harmonic frequencies exactly.'))
    content.append(f"<p>Across the six bands: RMSE {v['RMSE_cm1']:.2f} cm⁻¹, MAE {v['MAE_cm1']:.2f} cm⁻¹, mean absolute relative error {v['MAPE_percent']:.2f}%. The largest isotope redshift error is {100*physics['vibration_gate']['max_isotope_shift_fraction_error']:.3f} percentage points. All three advance criteria pass: maximum band-position error below 5%, shift error below 2 percentage points, and numerical spread below 1 cm⁻¹. The largest spread across three starting geometries and three Hessian displacements is {physics['vibration_gate']['max_numerical_spread_cm1']:.3f} cm⁻¹. These are pilot thresholds, not spectroscopic accuracy.</p>")
    content.append(img(root,'IR_band_positions','Band-position comparison only. Equal-height markers do not represent calculated or measured intensities. No full spectral line-shape comparison was performed.'))
    content.append(f"<p>The secondary GFN1-xTB calculation has a maximum error of {v['GFN1_max_error_percent']:.2f}%, failing the same 5% screen. The primary result was not replaced by the better-looking member of a method search. Adding HDO after primary validation, with the identical Hessian and no refitting, gives a maximum error of {extra['max_absolute_error_percent']:.2f}%. Its observed bend, OD-stretch and OH-stretch positions are 1402.20, 2726.73 and 3707.47 cm⁻¹. {ref('nist_hdo')}</p>")
    content.append(table(['HDO mode','Calculated / cm⁻¹','Observed / cm⁻¹','Error / %'],[[m,f'{p:.2f}',f'{o:.2f}',f'{e:+.3f}'] for m,p,o,e in zip(extra['mode_order'],extra['predicted_cm1'],extra['observed_IR_cm1'],extra['error_percent'])]))
    content.append(f"<p>Two execution failures were retained in the audit: the first optimization tolerance exceeded the 1 cm⁻¹ repeatability threshold, and this Windows xTB build produced identical built-in H/D frequencies despite reading the isotope input. The revised run uses tighter optimization and independently diagonalizes the saved Hessian. It reproduces xTB's default-mass frequencies within 0.005 cm⁻¹; the H and D electronic Hessians are identical. Multiplying all masses by four halves the frequencies exactly. H/D label exchange changes HDO results by only {extra['H_D_exchange_max_error_cm1']:.2g} cm⁻¹. A 101-point fictitious mass sweep is monotonic. See amendment.json and {ref('xtb_issue','upstream issue')}.</p>")
    content.append('<p><strong>What quantum mechanics adds:</strong> the electronic method supplies force constants from an approximate quantum electronic model; a quantum oscillator adds discrete transition energies and zero-point motion. However, ω² is an eigenvalue of M⁻¹ᐟ²KM⁻¹ᐟ² in both classical and quantum harmonic mechanics. Thus, the isotope frequency shift itself does not distinguish quantum nuclei from classical masses. Errors include approximate electronic forces and comparing harmonic frequencies with anharmonic fundamentals. Experimental tabulation uncertainty is about 0–1 cm⁻¹; model errors are tens to over 100 cm⁻¹. GFN methods are semiempirical: historical training overlap was not audited, so this is not a certified blind holdout.</p>')
    content.append(f"<p><strong>A stricter observable check exposes error cancellation.</strong> NIST also tabulates H₂O harmonic reference values of 1649, 3832 and 3943 cm⁻¹. Against those like-for-like harmonic references, the largest GFN2 error is {v['H2O_harmonic_reference_max_error_percent']:.2f}%, which exceeds 5%. The favorable fundamental-position comparison therefore does not imply equally accurate force constants: underestimated harmonic frequencies partially compensate for omitted anharmonic redshifts. This post-primary diagnostic limits the interpretation of the original pass. The optimized angle is about 107.22°, versus the listed equilibrium value 104.4776°. {ref('harmonic')}</p>")
    content.append('<h2>2. Liquid organization: a controlled exploratory test</h2><p>The established flexible q-TIP4P/F interaction model is held fixed while changing isotope masses and nuclear treatment. Twelve main runs cover H₂O/D₂O × classical/32-bead quantum × three independent orientation, velocity and thermostat seeds. Each has 64 molecules, 5 ps equilibration, 5 ps structural sampling and a 5 ps unthermostatted mobility branch. All use the same number density and a 0.5 fs integration step. The oxygen packing begins from one common template; these are independent stochastic trajectories, not independent preparation histories. '+ref('qtip','Model and parameters')+'; '+ref('rpmd','sampling method')+'.</p><p>“Beads” are imaginary-time replicas used to sample nuclear quantum equilibrium, not independent observations. Structure is averaged over beads; mobility uses centroid ring-polymer trajectories with the thermostat disabled. RPMD real-time motion is approximate. P=1 is the classical counterpart. The four-way isotope × nuclear-treatment design separates classical mass-dependent dynamics from a nuclear quantum contribution.</p>')
    content.append(table(['Condition','Bond mean / Å','H-bond degree','Tetrahedral order','Neighbors replaced at 1 ps'],[[name,ci(g['bond_mean_A']),ci(g['hydrogen_bond_degree'],3),ci(g['tetrahedral_order'],3),ci(g['neighbor_replaced_fraction_1ps'],3)] for name,g in groups.items()]))
    content.append('<p>Brackets are Student-t 95% intervals across three independent runs (two degrees of freedom); they do not include force-field error, finite size, equilibration bias or bead truncation. H-bond degree counts donor plus acceptor bonds using O–O &lt;3.5 Å and donor alignment within 30°. Neighbor replacement is a short-time structural diagnostic.</p>')
    content.append(img(root,'molecular_structure','Top-row intervals use independent runs. Neutron bands denote quoted ±0.005 Å uncertainties; their experimental estimator differs from the bare simulated radial peak. Temperature pilots have one seed each and cannot establish robust trends.'))
    contrast=md['paired_contrasts']
    content.append(table(['Contrast','Bond mean change / Å','H-bond degree change'],[[name,ci(contrast['bond_mean_A'][name]),ci(contrast['hydrogen_bond_degree'][name],3)] for name in ['classical_D_minus_H','quantum_D_minus_H','NQE_H','NQE_D','isotope_by_quantum_interaction']]))
    content.append('<p>At fixed potential, volume and temperature, the classical canonical coordinate distribution is proportional to exp(−V/kBT). Integrating momenta removes masses from the normalized coordinate distribution. A persistent equilibrium structure isotope effect therefore needs a nonclassical contribution or a changed potential/state; differences in these short classical samples can be sampling bias. Classical dynamical isotope effects remain possible. The interaction contrast is [(D−H)quantum − (D−H)classical], with pairing by seed.</p>')
    classical_fail=[key for key in ['bond_mean_A','hydrogen_bond_degree','tetrahedral_order'] if np.prod(contrast[key]['classical_D_minus_H']['CI95'])>0]
    if classical_fail:content.append('<p><strong>Equilibrium-null warning:</strong> the nominal 95% intervals exclude zero for classical D−H differences in '+', '.join(classical_fail)+'. Because exact classical equilibrium has no such mass dependence at this fixed state, this is evidence that short-run bias or sampling uncertainty is not adequately controlled. These are exploratory, unadjusted multiple comparisons. They strengthen the case against accepting a bulk quantum-improvement claim.</p>')
    content.append(table(['Condition','Bond-distribution width / Å','Mean-bond NQE change / Å'],[[name,ci(g['bond_distribution_sd_A']),ci(contrast['bond_mean_A']['NQE_H' if name.startswith('H2O') else 'NQE_D']) if name.endswith('P32') else 'Classical baseline'] for name,g in groups.items()]))
    content.append('<p>These widths describe pooled bond-length distributions, not uncertainty in their means. Quantum sampling can broaden the spatial distribution and change the mean on an anharmonic potential. Those are additions beyond changing classical masses; whether their predicted magnitude is accurate still needs observable-matched experimental validation.</p>')
    content.append('<p>An empirical classical potential can absorb quantum effects into fitted parameters. Holding the potential fixed here isolates explicit nuclear quantum treatment; it does not compare against every possible reparameterized classical model. Once density and viscosity are measured, the continuum flow calculations need no explicit nuclear quantum dynamics. The value of a microscopic quantum model would be accurate, transferable prediction of those properties before measuring them.</p>')
    peakrows=[]
    for name,g in groups.items():
        observed=.990 if name.startswith('H2O') else .985
        peakrows.append([name,f"{g['intramolecular_RDF_peak_A']['mean']:.4f}",f'{observed:.3f} ± 0.005',f"{g['intramolecular_RDF_peak_A']['mean']-observed:+.4f}"])
    content.append(table(['Condition','Simulated radial peak / Å','Neutron distance / Å','Difference / Å'],peakrows))
    content.append('<p>The neutron benchmark at 300 K gives approximately 0.990(5) Å for O–H and 0.985(5) Å for O–D. These are features of isotope-difference diffraction distributions, not mean bond lengths. We compare an intramolecular radial peak, smoothed with a specified two-bin Gaussian, and separately report means. We have not simulated the complete isotope-difference scattering signal or instrument truncation, so this is an approximate diagnostic, not a like-for-like statistical validation. The model was developed using ordinary-water measurements; “quantum improves agreement” is not automatically independent evidence. '+ref('neutron')+'.</p>')
    content.append(img(root,'molecular_stress','Failed checks remain visible. Neither a favorable short-time diffusion slope nor a local structural match is accepted as a validated transport coefficient.'))
    content.append(f"<p>Of the 12 main mobility runs, {mstats['main_window_passes']} pass the &lt;20% agreement screen between 1–2 ps and 2–4 ps diffusion fits, and {mstats['main_energy_passes']} pass the absolute energy-drift screen of 0.005 NkBT over 5 ps. The total ring energy is divided by bead count before this normalization. {mstats['structural_controls_passes']} of {mstats['structural_controls_total']} single-run structural sensitivity controls meet both specified tolerances (0.002 Å in mean bond length and 0.15 in H-bond degree). Passing an individual control does not establish convergence.</p>")
    content.append('<p>The predeclared energy screen uses final-minus-initial energy. This can include bounded integration oscillations and does not by itself identify secular drift. The saved summaries also give the energy range and a fitted time slope. The screen is retained unchanged rather than relaxed after observing failures.</p>')
    content.append(table(['Control','Bond change / Å','H-bond degree change','Combined screen'],[[r['tag'],f"{r['bond_difference_A']:+.5f}",f"{r['hbond_difference']:+.4f}",'PASS' if r['structural_screen_pass'] else 'FAIL'] for r in md['controls']]))
    small=[r for r in runs if r['dt_ps']==.00025]
    if small:
        s=small[0];content.append(f"<p>The added 0.25 fs trajectory gives an energy drift of {s['NVE_energy_drift_NkBT']:+.4f} NkBT ({'PASS' if s['NVE_energy_drift_pass'] else 'FAIL'}) and diffusion-window disagreement {100*s['diffusion_window_relative_disagreement']:.1f}% ({'PASS' if s['diffusion_window_pass'] else 'FAIL'}). It remains one trajectory, without a converged size-and-duration series for diffusion.</p>")
    content.append(f"<p>Controls also use 16/64 beads, 216 molecules and temperatures 280/320 K. Changing box size also changes the cutoff here, so that control combines size and cutoff sensitivity. All {len(runs)} runs completed {mstats['total_simulated_ps']:.0f} ps in total, counting equilibration and all branches once per trajectory, not once per bead. Potential energies and half-trajectory structural differences are available in raw summaries. No production viscosity or validated bulk diffusion value is claimed.</p>")
    content.append('<p>For context only, Mills reports self-diffusion values of 2.299 and 1.872 ×10⁻⁹ m²/s at 298.15 K for H₂O and D₂O. Our runs are at 300 K and fail transport-convergence requirements, so these values were not used to score an improvement. '+ref('mills')+'.</p>')
    content.append(f"<p>Implementation verification is separate from model validation: an independent pair-energy calculation agrees within {force['energy_error']:.2g} kJ/mol and finite-difference forces within {force['max_force_error']:.2g} kJ mol⁻¹ nm⁻¹. An exactly solvable harmonic oscillator checks the quantum sampler at P=1,16,32,64. Sample variance errors against finite-bead theory range from −3.34% to +2.71%; P=32 itself underestimates the exact continuous quantum variance by about 3.6% at 3700 cm⁻¹ and 300 K. The sampler check does not validate the water force field.</p>")
    content.append('<h2>3. Measured material properties to Navier–Stokes flow</h2><p>We use measurement-based IAPWS density and viscosity correlations at atmospheric pressure. These are external constitutive inputs; neither the xTB calculation nor the short molecular runs predict them. '+ref('iapws')+'.</p>')
    content.append(table(['25 °C property','H₂O','D₂O'],[['Density / kg m⁻³',f"{real['H2O']['rho']:.4f}",f"{real['D2O']['rho']:.4f}"],['Dynamic viscosity / mPa s',f"{1000*real['H2O']['mu']:.6f}",f"{1000*real['D2O']['mu']:.6f}"],['Steady mean channel speed / mm s⁻¹',f"{1000*real['H2O']['steady_mean_ms']:.6f}",f"{1000*real['D2O']['steady_mean_ms']:.6f}"],['Startup time constant / s',f"{real['H2O']['tau_s']:.6f}",f"{real['D2O']['tau_s']:.6f}"],['Dissipation after 1 s / W m⁻²',f"{real['H2O']['dissipation_Wm2']:.7g}",f"{real['D2O']['dissipation_Wm2']:.7g}"]]))
    content.append('<p>The matched channel has a 1 mm full gap, no-slip walls, initially stationary liquid and pressure gradient G=10 Pa/m for the table. The exact relations are ū∞=GH²/(12μ), τ=ρH²/(π²μ), and integrated steady dissipation per wall area =GHū. Viscosity determines steady speed; density changes acceleration and relaxation. At fixed pressure gradient D₂O dissipates less power because it flows more slowly. At fixed imposed speed the comparison would differ.</p>')
    content.append(f"<p>There are 36 main channel cases spanning 10/25/40 °C, G=1/10/100 Pa/m and four property combinations: ordinary water, D₂O density alone, D₂O viscosity alone, and both. The largest profile error against the analytic startup solution is {flowstats['channel_worst_relative_error']:.2g}; the largest discrete energy-balance residual is {flowstats['channel_worst_energy_residual']:.2g}. Convergence runs vary grid and time step separately. These are numerical verification against an exact solution, not new experimental velocity measurements.</p>")
    content.append(img(root,'fluid_connection','Matched continuum comparisons. The historical capillary comparison tests material-property consistency; it is not an independent velocity-field experiment.'))
    eb=read(root/'channel-error-budget.json')['comparisons']
    hb=[r for r in eb if r['isotope']=='H2O'];first,last=hb[0],hb[-1]
    content.append(f"<p>A separate error analysis explains why finer settings do not guarantee a monotonic reduction in total error. In the 129-node H₂O channel at 0.1 s, reducing the time step from 1 to 0.25 ms reduces the time-integration error from {100*first['temporal_error_relative_to_continuum']:.6f}% to {100*last['temporal_error_relative_to_continuum']:.6f}%, yet total error increases from {100*first['total_error_relative_to_continuum']:.6f}% to {100*last['total_error_relative_to_continuum']:.6f}%. Exact discrete-Laplacian modes show that spatial and temporal errors point in almost opposite directions and partly cancel. The smaller time step removes part of that cancellation while leaving spatial error unchanged. The reconstructed errors agree with the recorded solver results within {max(r['reconstructed_vs_recorded_error_difference'] for r in eb):.2g}.</p>")
    content.append(f"<p>The three-dimensional test reuses the project's incompressible Fourier solver unchanged. It evolves a periodic Taylor–Green vortex with the same initial velocity, domain and zero forcing for both liquids. The original 48 cases include three temperatures, three initial speed scales, property ablations, grid/time-step controls and three independent 1% perturbations for each isotope. Another 18 runs refine the strongest original speed across temperature. The maximum final divergence is {flowstats['vortex_worst_divergence']:.2g} in dimensionless units, and the maximum original energy residual is {flowstats['vortex_worst_energy_residual']:.2g}. An exact two-dimensional Taylor–Green solution provides a separate solver check.</p>")
    content.append('<p>The periodic cube has side 2π mm (approximately 6.283 mm); 1 mm is the length scale used to nondimensionalize its equations. The coarsest grid fails the 1% comparison screen in some stronger cases. Across the original sweep, the largest 24³-to-32³ field difference is '+f"{100*flowstats['strongest_forcing_worst_fine_field_difference']:.3f}%"+'. The highest original speed therefore uses 32³ fields in the isotope-difference summary. Stability tests measure finite-time perturbation amplification in this chosen family. Channel perturbations are restricted to unidirectional profiles; they do not establish general hydrodynamic stability, transition thresholds or global Navier–Stokes regularity.</p>')
    content.append(table(['T / K','Speed scale / mm s⁻¹','Grid','RMS field difference / mm s⁻¹','Difference relative to H₂O / %'],[[r['T_K'],r['U_scale_ms']*1000,r['grid_n'],f"{1000*r['rms_velocity_difference_ms']:.4f}",f"{100*r['relative_velocity_field_difference']:.3f}"] for r in flowstats['vortex_velocity_differences']]))
    gains=[r['absolute_perturbation_amplitude_gain'] for r in flowstats['perturbation_controls']]
    relative_gains=[r['relative_to_base_amplitude_gain'] for r in flowstats['perturbation_controls']]
    content.append(f"<p>For six perturbed 25 °C runs at 10 mm/s, final-to-initial absolute perturbation amplitude ratios range from {min(gains):.3f} to {max(gains):.3f}. Relative to the decaying base-flow amplitude, their gains range from {min(relative_gains):.3f} to {max(relative_gains):.3f}; a perturbation can shrink absolutely while growing as a fraction of the base flow. This is a finite-time result at 0.1 s, not a universal stability statement.</p>")
    content.append('<h2>Additional stress: push the flow beyond the original range</h2><p>After the initial results, the speed scale was raised to 50 and 100 mm/s, with 24³, 32³ and 48³ grids, both liquids at 25 °C, and a fixed 0.05 s duration. The advance screen remained a &lt;1% change in the full velocity field. These 12 runs deliberately bracket where the available resolution becomes inadequate.</p>')
    content.append('<p>The reported percentage is 100 × ‖u<sub>coarse</sub> − u<sub>fine</sub>‖₂ / ‖u<sub>fine</sub>‖₂ at the final time, using Fourier coefficients and including modes absent from the coarse grid. It is an overall field comparison, not a uniform percentage error at every point and not an isotope difference. The speed is the initial velocity scale; the vortex subsequently decays.</p>')
    content.append(table(['Liquid','Speed / mm s⁻¹','Grid pair','Field change / %','1% screen'],[[r['isotope'],1000*r['U_ms'],f"{r['coarse_n']}³ → {r['fine_n']}³",f"{100*r['relative_field_change']:.3f}",'PASS' if r['passed'] else 'FAIL'] for r in stress]))
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for U,style in [(.05,'-'),(.1,'--')]:
        for iso in COLORS:
            rows=[r for r in stress if r['U_ms']==U and r['isotope']==iso]
            axes[0].semilogy([r['fine_n'] for r in rows],[100*r['relative_field_change'] for r in rows],style+'o',color=COLORS[iso],label=f'{iso}, {1000*U:g} mm/s')
    axes[0].axhline(1,color='black',ls=':',label='1% screen');axes[0].set(xticks=[32,48],xticklabels=['24³ → 32³','32³ → 48³'],ylabel='Velocity-field difference (%)',title='Higher speed exposes insufficient resolution');axes[0].legend(fontsize=8)
    values=[v['max_relative_error_percent'],extra['max_absolute_error_percent'],v['GFN1_max_error_percent'],v['H2O_harmonic_reference_max_error_percent']]
    axes[1].bar(np.arange(4),values,color=['#176d9c','#176d9c','#cf6b2f','#cf6b2f']);axes[1].axhline(5,color='black',ls=':');axes[1].set(xticks=np.arange(4),xticklabels=['GFN2 H/D\nfundamentals','GFN2 HDO\nfundamentals','GFN1 H/D\nfundamentals','GFN2 H₂O\nharmonics'],ylabel='Largest frequency error (%)',title='The apparent accuracy depends on the test');axes[1].tick_params(axis='x',labelsize=8)
    savefig(fig,root,'stress_boundaries');content.append(img(root,'stress_boundaries','Both numerical resolution and the choice of physical observable matter. The harmonic-reference comparison is a post-primary diagnostic; it does not erase the original fundamental-position result.'))
    content.append('<p>A failed resolution test is a limit of this computation, not evidence of a physical singularity or a new quantum fluid effect. Small energy residuals alone do not establish spatial convergence.</p>')
    content.append('<h2>Follow-up: does the 100 mm/s result settle on finer grids?</h2><p>The failed 32³-to-48³ check motivated a bounded refinement to 64³ and 96³, plus a half time step at 64³ for each liquid. Geometry, initial flow, material properties and final physical time remain identical. The grid screen stays at 1%; the separate time-step screen is 0.1%. These six additional runs preserve the earlier failure.</p>')
    content.append(table(['Liquid','Comparison','Field change / %','Screen'],[[r['isotope'],f"{r['coarse_n']}³ → {r['fine_n']}³" if r['kind']=='grid' else '64³: Δt 0.01 → 0.005',f"{100*r['relative_field_change']:.6g}",'PASS' if r['passed'] else 'FAIL'] for r in refinement]))
    fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained')
    for iso in COLORS:
        rr=[r for r in stress if r['isotope']==iso and r['U_ms']==.1]+[r for r in refinement if r['isotope']==iso and r['kind']=='grid']
        ax.semilogy(range(len(rr)),[100*r['relative_field_change'] for r in rr],'o-',color=COLORS[iso],label=iso)
    ax.axhline(1,color='black',ls=':',label='1% comparison screen')
    ax.set(xticks=range(4),xticklabels=['24³ → 32³','32³ → 48³','48³ → 64³','64³ → 96³'],ylabel='Velocity-field difference (%)',title='100 mm/s: successive grid comparisons')
    ax.legend();savefig(fig,root,'flow_refinement')
    content.append(img(root,'flow_refinement','The ordinate compares two numerical fields for the same liquid. It is neither an H₂O–D₂O difference nor a measured error against the true flow.'))
    fh=np.load(root/'fluid-refinement/H2O-n96-dt0.01.npz')['final_hat'];fd=np.load(root/'fluid-refinement/D2O-n96-dt0.01.npz')['final_hat']
    fw=np.full((1,1,49),2.);fw[...,0]=fw[...,-1]=1.
    fnorm=lambda x:np.sqrt(np.sum(abs(x)**2*fw)/96**6)
    fast=dict(T_K=298.15,U_scale_ms=.1,final_time_s=.05,n=96,H2O_rms_velocity_ms=.1*fnorm(fh),D2O_rms_velocity_ms=.1*fnorm(fd),relative_material_field_difference=fnorm(fd-fh)/fnorm(fh),relative_RMS_speed_reduction=1-fnorm(fd)/fnorm(fh))
    ch=np.load(root/'fluid-refinement/H2O-n64-dt0.01.npz')['final_hat'];cd=np.load(root/'fluid-refinement/D2O-n64-dt0.01.npz')['final_hat'];cs=Solver(64,1.,workers=1)
    fast['material_field_difference_on_n64']=np.sqrt(cs.inner(cd-ch,cd-ch)/cs.inner(ch,ch))
    fast['isotope_difference_field_grid_sensitivity']=difference(cd-ch,fd-fh,cs)
    dump(root/'fast-flow-material-comparison.json',fast)
    content.append(f"<p>For comparison, using the same 96³ grid for both liquids gives an H₂O–D₂O velocity-field difference of {100*fast['relative_material_field_difference']:.3f}% relative to H₂O, with final RMS speeds {1000*fast['H2O_rms_velocity_ms']:.3f} and {1000*fast['D2O_rms_velocity_ms']:.3f} mm/s. That is the material-property comparison. It is distinct from the successive-grid differences above. In this unforced constant-density velocity equation, the material dependence enters through kinematic viscosity μ/ρ; these inputs still come from measured-property correlations.</p>")
    content.append(f"<p>The smaller isotope-difference field deserves its own precision check. Its magnitude changes only from {100*fast['material_field_difference_on_n64']:.4f}% to {100*fast['relative_material_field_difference']:.4f}% of the H₂O field between 64³ and 96³. However, the difference field itself changes by {100*fast['isotope_difference_field_grid_sensitivity']:.2f}% relative to its own smaller norm. Thus the overall-flow 1% screen passes, while a 1%-accurate spatial map of the isotope effect has not been established. This is an additional diagnostic of saved fields, not a replacement for the predeclared screen.</p>")
    content.append('<p>A finer grid remains an approximation. Two grids can agree while sharing a modeling error or missing the same physics. Shrinking differences support numerical convergence within these equations and this setup; they do not prove the finest field is exact, or substitute for an experimental velocity-field benchmark.</p>')
    content.append('<h2>Experimental comparison, uncertainty and provenance</h2>')
    content.append(table(['Temperature / °C','Measured μD/μH','IAPWS μD/μH','Error / %'],[[r['T_C'],r['observed_viscosity_ratio'],f"{r['predicted_viscosity_ratio']:.6f}",f"{r['relative_error_percent']:+.3f}"] for r in bench]))
    content.append(f"<p>Historical capillary measurements and the modern correlation agree within {flowstats['historical_max_viscosity_ratio_error_percent']:.3f}% at these three temperatures. However, those historical measurements are included in the correlation's literature/data provenance, so this is not holdout validation. The historical sample also contained excess oxygen-18. Inverse viscosity ratios give measurement-derived capillary conductance ratios, not separately measured velocity fields. The IAPWS implementation also reproduces all seven official heavy-water viscosity verification values within 10⁻⁷ relative error. {ref('hardy')}; {ref('visc_fit')}; {ref('visc')}.</p>")
    later=read(root/'later-viscosity-comparison.json');lr=later['comparisons']
    content.append(f"<p>A later check uses 2022 Brownian-motion viscosity measurements: all non-calibration liquid-state points above 277 K from H₂O run 6 and the two highest-purity heavy-water runs. The maximum discrepancy is {later['max_absolute_error_percent']:.3f}%. Heavy samples contain 96.66 mol% D₂O; predictions apply the paper's approximate linear mixture rule. Measurements share an older 293.15 K calibration, so this checks later temperature dependence rather than an independent absolute scale. {ref('later_visc')}</p>")
    content.append(table(['Run','T / K','Measured / mPa s','Predicted / mPa s','Error / %'],[[r['run'],r['T_K'],f"{r['observed_mean_mPas']:.4f}",f"{r['predicted_mPas']:.4f}",f"{r['relative_error_percent']:+.3f}"] for r in lr]))
    csvfile(root/'later_viscosity_comparison.csv',[{k:v for k,v in r.items() if k!='readings_mPas'} for r in lr])
    fig,ax=plt.subplots(figsize=(7,4.6),layout='constrained')
    for r in lr:
        ax.errorbar(r['observed_mean_mPas'],r['predicted_mPas'],xerr=r['approx_uncertainty_mPas'],fmt='o',color=COLORS['H2O' if r['x_D2O']==0 else 'D2O'],capsize=4)
    ax.plot([1.05,2.02],[1.05,2.02],color='gray',lw=1);ax.set(xlabel='Later measured viscosity (mPa s)',ylabel='IAPWS / mixture prediction (mPa s)',title='Later data agree; calibration is shared')
    savefig(fig,root,'later_viscosity');content.append(img(root,'later_viscosity','Blue: H₂O. Orange: 96.66% D₂O samples. Diagnostic bars combine 1.5% measurement uncertainty, assumed 0.15 K temperature uncertainty, and heavy-sample composition uncertainty of 0.04. Shared calibration prevents treating the points as independent absolute tests.'))
    q=sens['steady_velocity_D_over_H_quantiles'];tq=sens['decay_time_D_over_H_quantiles']
    content.append(f"<p>An explicitly assumed sensitivity calculation varies H₂O viscosity ±0.5%, D₂O viscosity ±1%, and both densities ±0.02%, independently and uniformly, for 100,000 draws. Its central 95% range is [{q[0]:.5f}, {q[2]:.5f}] for steady D/H speed and [{tq[0]:.5f}, {tq[2]:.5f}] for D/H relaxation time. These are sensitivity ranges, not calibrated confidence intervals for the experiments or IAPWS correlations.</p>")
    content.append('<p>A further literature check found H₂O/D₂O capillary experiments with silica colloids and velocity profiles inferred from chromatography peaks (Yoon &amp; Yoon, 2018). The authors interpret temperature anomalies as quantum behavior at interfaces; they report no corresponding anomaly without silica. That setup includes colloids, interfaces, tracer transport and imposed pump flow, which this pure-fluid bulk model does not represent. It was therefore not scored as validation of our pressure-driven channel or periodic vortex, nor counted as proof that our simulated flow difference requires quantum nuclei. '+ref('interface_flow')+'.</p>')
    content.append('<p>Existing project work informed the design: the isotope-mass pilot separates short-time inertial response from structure; the molecule-influence experiment motivates neighbor-replacement diagnostics; the fluid-organization pilot cautions against interpreting passive markers as molecules. Only the Navier–Stokes solver is reused directly. It is pinned to repository snapshot 4e4a970e20d3289e883dec5c34e1a1bd3968795e; SHA-256 is recorded in provenance.json. Prior reported results were not rerun or counted among this study\'s completed tests. '+ref('repo')+'.</p>')
    content.append('<h2>What remains unestablished</h2><ul><li>A quantitatively validated liquid nuclear-quantum improvement across temperature and preparation.</li><li>Converged diffusion and viscosity from the same molecular model, including finite-size, long-time, bead and timestep checks.</li><li>A direct quantum-calculation → viscosity → flow prediction. This pilot links independently measured material properties to flow.</li><li>An independent experimental velocity-field benchmark with matched H₂O/D₂O geometry and forcing.</li></ul><p>The next expansion should first lengthen and refine the failing molecular controls, forward-model the neutron observable and obtain converged transport coefficients. Only then should any model-predicted viscosity be inserted into the fluid solver and judged against an independent flow experiment. The completed package is an executable, falsifiable pilot with explicit failures, not a claim that all three research goals are validated.</p>')
    content.append('<h2>Reproduce and inspect</h2><p>README.txt gives installation, verification and full rerun commands. run_research.py refuses to overwrite an existing results directory. Raw Hessians, logs, molecular trajectories, saved fluid fields, protocols, amendment records, numerical tables and PNG/SVG charts are included. verify_results.py checks checksums and scientific invariants without rerunning the expensive simulations. GPU stochastic trajectories can differ across platforms; statistical agreement is the reproduction target, not bitwise identity.</p>')
    content.append('<h2>Primary sources</h2><ol>'+''.join(f'<li id="{id}"><a href="{url}">{html.escape(title)}</a></li>' for id,title,url in SOURCES)+'</ol>')
    css='''body{margin:0;background:#edf1f2;color:#243746;font:17px/1.6 system-ui,Segoe UI,sans-serif}main{max-width:1100px;margin:0 auto;background:white;padding:56px 60px}header{border-bottom:3px solid #176d9c;padding-bottom:30px}.eyebrow{font-size:12px;letter-spacing:1.3px;color:#507183}h1{font-size:46px;line-height:1.13;letter-spacing:-1.3px;margin:18px 0}h2{margin-top:48px;font-size:28px;line-height:1.25}p{max-width:98ch}.lede{font-size:21px;color:#486372}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:32px 0}.cards div{background:#eff5f7;padding:22px;border-top:3px solid #176d9c}.cards b{display:block;font-size:32px}.cards span{font-size:14px;display:block}.table{overflow-x:auto;margin:24px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:11px 12px;border-bottom:1px solid #dce3e7;vertical-align:top}th{background:#edf3f5}tr:nth-child(even){background:#f9fbfc}figure{margin:32px -24px}img{width:100%;height:auto}figcaption{font-size:13px;color:#586b76;padding:5px 24px}a{color:#176d9c}li{margin:8px 0}code{font-size:14px}@media(max-width:700px){main{padding:24px}h1{font-size:34px}.cards{grid-template-columns:1fr}figure{margin:20px -10px}}@media print{body{background:white}main{padding:0}h2{break-after:avoid}figure,table{break-inside:avoid}.cards{grid-template-columns:repeat(3,1fr)}}'''
    (root/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Quantum water: completed tests and stress limits</title><style>'+css+'</style><main>'+''.join(content)+'</main></html>',encoding='utf-8')
    dump(root/'sources.json',[dict(id=i,title=t,url=u) for i,t,u in SOURCES])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);a=ap.parse_args()
    vr,sh,vs=vibration(a.data);fs,re,be=flow(a.data);md,ms=molecular(a.data)
    report(a.data,vr,sh,vs,fs,re,be,md,ms)
    print('Built report.html, numerical tables and scientific figures.',flush=True)
