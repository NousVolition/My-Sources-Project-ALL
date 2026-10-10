"""Build the reviewed, self-contained scientific report and checksummed release."""
from pathlib import Path
from html import escape
import argparse,base64,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from water_battery import cases,save

HERE=Path(__file__).resolve().parent;D=HERE/'data'
COL={'water':'#087e8b','rotors':'#d27632'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.titleweight':'bold','svg.hashsalt':'water-battery-20261010'})


def fmt(value,digits=3):return '—' if value is None else f'{value:.{digits}f}'
def ms(obj):return fmt(obj['mean'])+(' ± '+fmt(obj['SD']) if obj['SD'] is not None else '')
def interval(obj):return '['+', '.join(fmt(x) for x in obj['nominal_95pct_t_interval'])+']' if obj['nominal_95pct_t_interval'] is not None else 'unavailable'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args()
    result=json.loads((D/('partial_results.json' if args.preview else 'results.json')).read_text())
    models=json.loads((D/'model_results.json').read_text());protocol=json.loads((D/'protocol.json').read_text());validation=json.loads((D/'integration_validation.json').read_text())
    if not args.preview:assert result['complete'] and len(result['refinements'])==12 and models['checks']['passed'] and validation['passed']
    target=HERE if not args.preview else HERE.parent.parent/'work'/'battery-preview';target.mkdir(exist_ok=True,parents=True)
    rows=result['rows'];aggregate={(r['kind'],r['case']):r for r in result['aggregates']};data={}
    for row in rows:
        key=(row['kind'],row['replica'],row['case']);data[key]=np.load(D/f'{key[0]}_r{key[1]}_{key[2]}.npz')
    figures=[]
    def finish(fig,name,title):
        for panel in fig.axes:
            if not panel.has_data():
                panel.text(.5,.5,'Results still running' if args.preview else 'No result available',ha='center',va='center',transform=panel.transAxes,color='#55707a')
        fig.suptitle(title,fontsize=16,x=.04,ha='left');fig.tight_layout(rect=(0,0,1,.96))
        fig.savefig(target/(name+'.png'),dpi=170,bbox_inches='tight');fig.savefig(target/(name+'.svg'),bbox_inches='tight');plt.close(fig);figures.append((name,title))
    def series(ax,name,key,kind='water',component=None,label=None,shift=0):
        ds=[data[kind,r,name] for r in range(3) if (kind,r,name) in data]
        if not ds:return
        y=np.array([d[key] if component is None else d[key][:,component] for d in ds]);t=ds[0]['time_ps']-shift
        color=COL[kind]
        for q in y:ax.plot(t,q,color=color,alpha=.2,lw=.6)
        ax.plot(t,y.mean(axis=0),color=color,lw=1.6,label=label or kind)

    fig,ax=plt.subplots(2,2,figsize=(12,8))
    for kind in ['water','rotors']:
        for amp,style in [(.1,'--'),(.5,'-')]:
            freq=[];response=[];sd=[];projection=[];psd=[]
            for f in [.1,.5,1.]:
                subset=[r for r in rows if r['kind']==kind and r['case']==f'rotate_a{amp:g}_f{f:g}']
                if not subset:continue
                values=[r['phase']['response_magnitude'] for r in subset];pp=[r['late_field_projection'] for r in subset]
                freq.append(f);response.append(np.mean(values));sd.append(np.std(values,ddof=1) if len(values)>1 else 0);projection.append(np.mean(pp));psd.append(np.std(pp,ddof=1) if len(pp)>1 else 0)
            label=f'{kind}, {amp:g} V/nm'
            ax[0,0].errorbar(freq,response,yerr=sd,color=COL[kind],ls=style,marker='o',capsize=3,label=label)
            ax[0,1].errorbar(freq,projection,yerr=psd,color=COL[kind],ls=style,marker='o',capsize=3)
    for a,title in zip(ax[0],['Phase-coherent orientation amplitude','Alignment with the instantaneous drive']):
        a.set_xscale('log');a.set(xlabel='Drive frequency (cycles/ps)',ylabel='Orientation (dimensionless)',title=title);a.grid(alpha=.2)
    ax[0,0].legend(fontsize=8)
    for j,kind in enumerate(['water','rotors']):
        matrix=np.full((2,3),np.nan)
        for i,amp in enumerate([.1,.5]):
            for k,f in enumerate([.1,.5,1.]):
                sub=[r for r in rows if r['kind']==kind and r['case']==f'rotate_a{amp:g}_f{f:g}']
                if sub:
                    count=sum(r['phase']['finite_window_tracking'] for r in sub);matrix[i,k]=count/len(sub);ax[1,j].text(k,i,f'{count}/{len(sub)}',ha='center',va='center',color='white' if count/len(sub)>.5 else '#102d3d',fontsize=18)
        ax[1,j].imshow(matrix,vmin=0,vmax=1,cmap='YlGnBu',alpha=.65,aspect='auto');ax[1,j].set(xticks=range(3),xticklabels=['0.1','0.5','1.0'],yticks=[0,1],yticklabels=['0.1','0.5'],xlabel='Drive frequency (cycles/ps)',ylabel='Field strength (V/nm)',title=kind+': runs meeting the tracking screen')
    finish(fig,'rotation_results','Following a rotating drive: response and loss of coordination')

    fig,ax=plt.subplots(1,3,figsize=(14,4.5))
    for j,name in enumerate(['ramp_fast','ramp_slow']):
        for kind in ['water','rotors']:
            ds=[data[kind,r,name] for r in range(3) if (kind,r,name) in data]
            for d in ds:ax[j].plot(d['field_V_nm'][:,0],d['polarization'][:,0],color=COL[kind],alpha=.2,lw=.7)
            if ds:ax[j].plot(ds[0]['field_V_nm'][:,0],np.mean([d['polarization'][:,0] for d in ds],axis=0),color=COL[kind],label=kind)
        ax[j].set(xlabel='Field along x (V/nm)',ylabel='Mean dipole direction along x',title=name.replace('_',' '));ax[j].grid(alpha=.2)
    for name,color in [('history_positive','#087e8b'),('history_negative','#963e7e')]:
        ds=[data['water',r,name] for r in range(3) if ('water',r,name) in data]
        if ds:
            for d in ds:ax[2].plot(d['time_ps']-10,d['polarization'][:,0],color=color,alpha=.2,lw=.7)
            ax[2].plot(ds[0]['time_ps']-10,np.mean([d['polarization'][:,0] for d in ds],axis=0),color=color,label=name.replace('history_','')+' preparation')
    ax[0].legend();ax[2].axvline(0,color='black',ls=':',lw=1);ax[2].set(xlabel='Time after common +0.1 V/nm field begins (ps)',ylabel='Mean dipole direction along x',title='Different histories, same final conditions');ax[2].legend(fontsize=8)
    finish(fig,'history_results','Finite-rate loops and the test for persistent history dependence')

    fig,ax=plt.subplots(2,2,figsize=(12,8));colors=['#b84237','#449045','#414bb0']
    for a,name in zip(ax.flat,['competing','chirp','sequence_xyz','sequence_scrambled']):
        ds=[data['water',r,name] for r in range(3) if ('water',r,name) in data]
        if ds:
            t=ds[0]['time_ps'];p=np.mean([d['polarization'] for d in ds],axis=0)
            for i,c in enumerate(colors):a.plot(t,p[:,i],color=c,label='xyz'[i]);a.plot(t,ds[0]['field_V_nm'][:,i],color=c,ls=':',alpha=.35,lw=.7)
        if name.startswith('sequence'):a.axvline(30,color='black',ls='--');a.text(31,.65,'drive off',fontsize=9)
        a.set(title=name.replace('_',' '),xlabel='Time (ps)',ylabel='Orientation / field (V/nm, dotted)',ylim=(-.8,.8));a.legend(fontsize=8,ncol=3)
    finish(fig,'complex_drives','Competing rhythms, changing frequency, and sequence training')

    fig,ax=plt.subplots(2,2,figsize=(12,8))
    for name,style in [('baseline',':'),('feedback_align','-'),('feedback_turn','--')]:
        ds=[data['water',r,name] for r in range(3) if ('water',r,name) in data]
        if ds:ax[0,0].plot(ds[0]['time_ps'],np.mean([np.linalg.norm(d['polarization'],axis=1) for d in ds],axis=0),ls=style,label=name.replace('_',' '))
    ax[0,0].set(title='Feedback based on collective orientation',xlabel='Time (ps)',ylabel='Collective orientation magnitude');ax[0,0].legend(fontsize=8)
    for rep in range(3):
        d=data.get(('water',rep,'feedback_turn'))
        if d is not None:ax[0,1].plot(d['polarization'][:,0],d['polarization'][:,1],alpha=.7,label=f'start {rep+1}')
    ax[0,1].set(title='Turning feedback: orientation trajectory',xlabel='Mean x direction',ylabel='Mean y direction');ax[0,1].set_aspect('equal',adjustable='datalim')
    if ax[0,1].lines:ax[0,1].legend(fontsize=8)
    for kind in ['water','rotors']:series(ax[1,0],'spatial_release','polarization',kind,0)
    ax[1,0].axvline(10,color='black',ls=':');ax[1,0].set(title='Spatial field removed at 10 ps',xlabel='Time (ps)',ylabel='Global x orientation');ax[1,0].legend(fontsize=8)
    series(ax[1,1],'thermal_quench','temperature_K');ax[1,1].axvline(10,color='black',ls=':');ax[1,1].set(title='Temperature pulse, then return to 300 K',xlabel='Time (ps)',ylabel='Temperature (K)')
    finish(fig,'engineered_environments','Unconventional environments: feedback, spatial structure, and heat')

    fig,ax=plt.subplots(1,3,figsize=(15,11));names=[c['name'] for c in cases()];y=np.arange(len(names))
    for kind,offset in [('water',-.12),('rotors',.12)]:
        for j,key in enumerate(['late_norm','late_spatial_excess']):
            vals=[aggregate.get((kind,n),{}).get(key,{}) for n in names]
            means=[v.get('mean',np.nan) for v in vals];errors=[v.get('SD') or 0 for v in vals]
            ax[j].errorbar(means,y+offset,xerr=errors,fmt='o',ms=3,capsize=2,color=COL[kind],label=kind)
    for i,name in enumerate(names):
        changes=[]
        for rep in range(3):
            a=next((r for r in rows if r['kind']=='water' and r['replica']==rep and r['case']==name),None)
            b=next((r for r in rows if r['kind']=='water' and r['replica']==rep and r['case']=='baseline'),None)
            if a and b:changes.append(a['late_hbond_neighbors']-b['late_hbond_neighbors'])
        if changes:ax[2].errorbar(np.mean(changes),i,xerr=np.std(changes,ddof=1) if len(changes)>1 else 0,fmt='o',ms=3,capsize=2,color=COL['water'])
    for j,a in enumerate(ax):
        a.set(yticks=y,yticklabels=[n.replace('_',' ') for n in names] if j==0 else [],ylim=(len(names)-.4,-.6));a.grid(axis='y',alpha=.17);a.axvline(0,color='#aaa',lw=.7)
    ax[0].set(title='Collective orientation magnitude',xlabel='Late mean |P|');ax[0].legend(fontsize=9)
    ax[1].set(title='Additional neighbor alignment',xlabel='Excess over shuffled-direction expectation')
    ax[2].set(title='Hydrogen-bond neighborhood change',xlabel='Late change versus undriven water')
    finish(fig,'all_protocols','All protocols: late observables, including small and absent effects')

    fig,ax=plt.subplots(1,3,figsize=(14,4.6));rec=result.get('recurrence')
    if rec:
        item=next(c for c in rec['clusters'] if c['k']==3);im=ax[0].imshow(item['transition'],vmin=0,vmax=1,cmap='Blues');fig.colorbar(im,ax=ax[0],shrink=.7)
        ax[0].set(title='Fitted 1 ps transition probabilities',xlabel='Next coarse state',ylabel='Current coarse state',xticks=range(3),yticks=range(3))
        for i,rep in enumerate([1,2]):
            vals=[r for r in rec['linear_forecasts'] if r['replica']==rep]
            for j,key in enumerate(['linear_prediction_MSE','persistence_MSE','constant_training_mean_MSE']):ax[1].plot([r['horizon_ps'] for r in vals],[r[key] for r in vals],marker='o',ls=['-','--',':'][j],color=['#087e8b','#963e7e'][i],label=f'start {rep+1}: '+['fitted','persistence','constant'][j])
        ax[1].set(title='Held-out orientation prediction',xlabel='Forecast horizon (ps)',ylabel='Mean squared error');ax[1].legend(fontsize=7)
        for i,r in enumerate(rec['GLV_forecasts']):
            con=rec.get('constrained_activity_forecasts',[{},{}])[i]
            vals=[r['GLV_MSE'],con.get('MSE'),r['persistence_MSE'],r['constant_training_mean_MSE']]
            for j,v in enumerate(vals):
                if v is not None:ax[2].bar(i*5+j,v,color=['#087e8b','#963e7e','#d27632','#818b91'][j])
        ax[2].set(title='Trying activity models on water data',ylabel='1 ps activity prediction MSE',xticks=[1.5,6.5],xticklabels=['held-out start 2','held-out start 3']);ax[2].text(.02,.98,'teal: unrestricted GLV\npurple: competition + source\norange: persistence; gray: constant',transform=ax[2].transAxes,va='top',fontsize=8)
    finish(fig,'reduced_models','Does a simpler description predict water trajectories?')

    ma=np.load(D/'model_trajectories.npz');fig,ax=plt.subplots(2,3,figsize=(14,8))
    ad=models['adler'];ax[0,0].plot([r['mu'] for r in ad],[r['analytic_mean_rate'] for r in ad],'-',label='analytic');ax[0,0].scatter([r['mu'] for r in ad],[r['late_regression_rate'] for r in ad],s=18,label='simulation');ax[0,0].set(title='Phase locking to phase drift',xlabel='Mismatch / coupling',ylabel='Mean phase rate');ax[0,0].legend(fontsize=8)
    for initial in ['rest','running']:
        sub=[r for r in models['josephson'] if r['beta']==10 and r['initial']==initial];ax[0,1].plot([r['current'] for r in sub],[r['mean_phase_rate'] for r in sub],'-o',label=initial)
    ax[0,1].set(title='Josephson model: history dependence',xlabel='Normalized current',ylabel='Mean phase rate');ax[0,1].legend(fontsize=8)
    for row in models['protein_feedback']:
        if row['n']==2:
            for root in row['roots']:ax[0,2].plot(row['beta'],root['p'],'.',ms=3,color='#087e8b' if root['stable'] else '#d27632')
    ax[0,2].set(title='Cooperative protein feedback',xlabel='Feedback strength',ylabel='Steady activity');ax[0,2].text(.03,.97,'teal: stable; orange: unstable',transform=ax[0,2].transAxes,va='top',fontsize=8)
    for eps in [.05,.1,.2,.5,1.]:
        z=ma[f'quench_{eps:g}'];ax[1,0].plot(z[:,0],z[:,1],label=f'{eps:g}')
    ax[1,0].set(title='Threshold delay from a small seed',xlabel='Time (model units)',ylabel='Amplitude',xlim=(0,230));ax[1,0].legend(fontsize=7,title='growth')
    for lam in [-.5,.1,1.]:
        z=ma[f'laser_{lam:g}_20'];ax[1,1].plot(z[:,0],z[:,1],label=f'pump {lam:g}');ax[1,1].plot(z[:,0],z[:,-1],ls=':',color='gray',lw=1)
    ax[1,1].set(title='Laser model and reduced approximation',xlabel='Time (model units)',ylabel='Field amplitude');ax[1,1].legend(fontsize=7)
    cusp=ma['cusp_stability'];im=ax[1,2].scatter(cusp[:,1],cusp[:,0],c=cusp[:,3],cmap='YlGnBu',s=8,vmin=0,vmax=2);ax[1,2].set(title='Cusp normal form: stable equilibria',xlabel='r',ylabel='k');fig.colorbar(im,ax=ax[1,2],shrink=.7,ticks=[0,1,2])
    finish(fig,'reference_thresholds','Mathematical reference tests — not molecular water models')

    fig,ax=plt.subplots(2,3,figsize=(14,8))
    for i,eta in enumerate([0.,.01]):
        z=ma[f'heteroclinic_{eta:g}'];ax[0,i].plot(z[:,0],z[:,1:]);ax[0,i].set(title='Three-state competition, noise '+str(eta),xlabel='Time (model units)',ylabel='Activity')
    z=ma['hierarchy_0.0001'];ax[0,2].plot(z[:,0],z[:,-3:]);ax[0,2].set(title='Nine-state model: imposed groups',xlabel='Time (model units)',ylabel='Group activity')
    z=ma['chain_pacemaker_1'];activity=z[:,1:].reshape(-1,16,3);ax[1,0].imshow(activity.argmax(axis=2).T,origin='lower',aspect='auto',extent=[0,z[-1,0],0,16],cmap='viridis',interpolation='nearest');ax[1,0].set(title='Directed chain with assigned first unit',xlabel='Time (model units)',ylabel='Unit')
    for beta in [.5,1.,2.]:
        z=ma[f'sir_{beta:g}'];ax[1,1].plot(z[:,0],z[:,2],label=f'contact {beta:g}')
    ax[1,1].set(title='Explicit SIR threshold reference',xlabel='Time (model units)',ylabel='Infected fraction');ax[1,1].legend(fontsize=8)
    for k in [-1.,0.,.5,2.]:
        z=ma[f'phase_all_{k:g}'];ax[1,2].plot(z[:,0],z[:,-2],label=f'coupling {k:g}')
    ax[1,2].set(title='Identical oscillators, different coupling',xlabel='Time (model units)',ylabel='Phase coherence');ax[1,2].legend(fontsize=8)
    finish(fig,'reference_networks','Switching, noise, coupling, and coarse descriptions')

    fig,ax=plt.subplots(1,2,figsize=(11,4.5))
    if 'attraction_vs_stability' in models:
        for row in models['attraction_vs_stability']['rows']:
            z=ma[f'attraction_{row["initial_angle"]:g}'];ax[0].plot(z[:,0],z[:,-1],label=f'initial angle {row["initial_angle"]:g}');ax[1].plot(z[:,0]/row['time_to_opposite_point'],z[:,-1],label=f'{row["initial_angle"]:g}')
        ax[0].set(title='Each trajectory eventually returns toward the point',xlabel='Time (model units)',ylabel='Distance from equilibrium on the circle');ax[0].legend(fontsize=8)
        ax[1].set(title='A small start can still make a large excursion',xlabel='Time / time to opposite point',ylabel='Distance from equilibrium on the circle');ax[1].axhline(1,color='gray',ls=':');ax[1].set_xlim(0,4)
    finish(fig,'attraction_reference','Attraction and stability are different measurements')

    def get(kind,name,key):return aggregate.get((kind,name),{}).get(key,dict(mean=None,SD=None))
    n_complete=result['completed_main_runs'];contrasts=result['contrasts'];history=contrasts.get('water_history_positive_minus_negative');slow=contrasts.get('water_slow_minus_fast_loop')
    positives=sum(r.get('phase',{}).get('finite_window_tracking',False) for r in rows if r['kind']=='water' and r['case'].startswith('rotate_'))
    rotorpositives=sum(r.get('phase',{}).get('finite_window_tracking',False) for r in rows if r['kind']=='rotors' and r['case'].startswith('rotate_'))
    statements=[f'{n_complete} of 156 planned main trajectories are complete. The full design contains 6,204 ps (6.204 ns) of molecular and noninteracting-control dynamics, plus 12 refinement trajectories and short numerical audits.',
        f'The six uniform rotating-field settings produce {positives} water and {rotorpositives} noninteracting-control runs meeting the finite-window tracking screen (18 runs per kind when complete). Low collective amplitude makes phase unreliable and is explicitly counted as failing this screen, rather than unwrapped through noise.',
        f'With a steady +0.5 V/nm field along x, late alignment is {ms(get("water","steady","late_field_projection"))} in water and {ms(get("rotors","steady","late_field_projection"))} without intermolecular interactions. Rotating-drive responses can differ in the opposite direction; the figures and every condition are retained.']
    rotation_comparison=[]
    for frequency in [.1,.5,1.]:
        name=f'rotate_a0.5_f{frequency:g}'
        water=[r for r in rows if r['kind']=='water' and r['case']==name]
        control=[r for r in rows if r['kind']=='rotors' and r['case']==name]
        if not water or not control:continue
        values=[f'{frequency:g}']
        for group in [water,control]:
            values += [fmt(np.mean([r['phase']['complex_response_real'] for r in group])),fmt(np.mean([r['phase']['response_magnitude'] for r in group])),fmt(np.rad2deg(np.mean([r['phase']['mean_lag_radians'] for r in group])),1)]
        rotation_comparison.append(values)
    if len(rotation_comparison)==3:
        slow_rows=[r for r in rows if r['kind']=='water' and r['case']=='rotate_a0.5_f0.1']
        lags=np.rad2deg([r['phase']['mean_lag_radians'] for r in slow_rows])
        statements.append(f'The steady-versus-rotating contrast holds in all three matched configurations at 0.5 V/nm. At the slowest rotation, water has a larger coherent response than the noninteracting control ({rotation_comparison[0][2]} versus {rotation_comparison[0][5]}), but lags the field by {lags.min():.1f}–{lags.max():.1f} degrees. Thus lower instantaneous alignment does not always mean weaker collective order. At the two faster rates the coherent water response also falls. This supports delayed collective response within this model; removing all interactions is not a sweep of coupling strength.')
        low=[r['late_temperature_K'] for r in rows if r['kind']=='water' and r['case']=='low_drag']
        if low:statements.append(f'The low-drag rotating runs reach {min(low):.0f}–{max(low):.0f} K despite a 300 K bath target. Drag changes both dissipation and temperature here, so those runs cannot isolate the effect of friction at a common temperature.')
    if history:statements.append(f'After opposite preparation histories and 50 ps under the same +0.1 V/nm field, the late water alignment difference is {ms(history)}, with nominal 95% t interval {interval(history)}. This finite-duration, three-start comparison does not establish equilibrium bistability.')
    if slow:statements.append(f'The slow-minus-fast absolute water loop area is {ms(slow)} (orientation × V/nm). Rate-dependent loops can arise from ordinary relaxation; a loop by itself is not proof of two stable water states.')
    sequence=[]
    for name in ['sequence_xyz','sequence_scrambled','baseline']:
        sub=[r['six_ps_coherence'] for r in result['sequence_recall'] if r['kind']=='water' and r['case']==name]
        if sub:sequence.append(name.replace('_',' ')+': '+fmt(np.mean(sub)))
    if sequence:statements.append('Teacher-off six-ps directional coherence (normalized Fourier amplitude): '+ '; '.join(sequence)+'. Correct tracking during training is not teacher-off recall; these are explicit control comparisons, not a learning claim.')
    if ('water','feedback_align') in aggregate and ('water','feedback_turn') in aggregate:
        statements.append('Alignment feedback gives late water |P|='+ms(get('water','feedback_align','late_norm'))+', versus '+ms(get('water','feedback_turn','late_norm'))+' for turning feedback and '+ms(get('water','baseline','late_norm'))+' without a drive. The alignment controller specifies no preferred spatial direction, but it deliberately supplies global feedback. The result belongs to that engineered water–controller system.')
    if result['refinements']:
        agreed=sum(r['phase']['finite_window_tracking']==r['base_tracking'] for r in result['refinements'])
        statements.append(f'The refined runs retain the original tracking decision in {agreed} of {len(result["refinements"])} comparisons. Response magnitudes and kinetic temperatures are recorded separately; stochastic trajectories are not expected to coincide point by point.')
    if rec:
        statements.append(f'Three principal components explain {100*rec["first_three_PC_variance_fraction"]:.1f}% of the training variance in 24 local orientation features. Coarse states are fitted on one trajectory and assessed on two held-out trajectories; the forecast plot includes persistence and constant baselines. The GLV fit is a deliberately speculative, axis-dependent mapping Aᵢ=Pᵢ², and poor predictions are retained.')
        ratios=[r['GLV_MSE']/r['persistence_MSE'] for r in rec['GLV_forecasts'] if r['GLV_MSE'] is not None]
        statements.append('The attempted competition-model forecast has '+', '.join(fmt(q,2)+'×' for q in ratios)+' the error of simply holding the last observed activity constant on the two held-out water trajectories. Ratios above one mean that this particular model attempt did worse than the simple baseline.')
        if rec.get('unrestricted_GLV_has_negative_competition_coefficients'):statements.append('The unrestricted GLV fit also has negative interaction coefficients, so it cannot be interpreted as pure competition. A separate fit enforces nonnegative competition and adds a nonnegative constant source; both comparisons remain exploratory, and both are reported.')
    if all('late_fraction_with_one_chunk_above_80pct' in r for r in models['engineered_hierarchy']['rows']):
        fractions=[r['late_fraction_with_one_chunk_above_80pct'] for r in models['engineered_hierarchy']['rows']]
        statements.append('The chosen nine-state modular reference model also failed to reproduce pronounced chunks: the late fractions with one group holding more than 80% of activity are '+', '.join(fmt(q) for q in fractions)+'. Its weaker group fluctuations are retained rather than retuned until they resemble the source figure.')
    method=[('What was changed','All 216 TIP3P molecules retain identical masses, charges, rigid geometry and pair interactions. The environment changes through uniform or spatially patterned fields, field histories, temperature or bath friction. The two feedback cases add an external controller based on collective orientation; they are engineered dynamical systems, not newly discovered intrinsic molecular laws.'),
        ('Controls and uncertainty','Each condition uses the same three saved starting configurations and separately seeded Langevin trajectories. Removing intermolecular forces gives overlapping rigid ghost dipoles, not a second water model. Frames and molecules are correlated; means and SD are across three run summaries. The intervals are exploratory, unadjusted for this broad search, and do not treat hundreds of frames as independent evidence.'),
        ('Driving can also heat the model','The bath target is 300 K except during the explicit 350 K pulse. Under continuous forcing, measured kinetic temperatures can exceed the bath target. The table reports them rather than assuming perfect thermal matching. Interaction changes can therefore alter both motion and heating; the present comparisons do not isolate these effects with a temperature-matched ensemble.'),
        ('Local structure','Additional neighbor alignment subtracts the exact expected dipole dot product after randomly reassigning the observed directions to the fixed oxygen positions: (N|P|²−1)/(N−1). It preserves the entire global orientation distribution. The resulting excess is a spatial association, not proof of causal cooperation. Neighbors use O–O<0.35 nm; hydrogen bonds also require a donor O–H angle within 30° of the donor-to-acceptor direction. Error bars in the all-protocols plot are between-run SD; all quantities there use the last third of each run, which may be a recovery window.'),
        ('Phase and waveform measurement','The rotating-field analysis begins after the larger of one period and one-quarter of the run. It requires |Pxy|≥0.15 on at least 95% of frames, a sufficiently long contiguous valid segment, phase concentration≥0.9, relative drift≤0.05 and phase excursion<π. Thresholds 0.1 and 0.2 are also reported. A failed screen can mean weak response, noisy orientation or phase slipping; it does not prove absence of every form of synchronization. Crossing this analysis cutoff is not a demonstrated saddle-node bifurcation: response amplitude may decrease smoothly with frequency. Field arrays report the nominal target at observation times; the implemented drive is its documented midpoint hold.'),
        ('History, sequence, and memory','Opposite preparations lead to the same final field for 50 ps. Fast and slow field loops use the same field range. Repeating x→y→z training is compared with an equally long scrambled sequence and a field-free control. External forcing is removed before recall measurement. The model has no changing molecular parameters or learned synaptic weights; finite residues, transients and imposed sequences are not permanent memory.'),
        ('Spatial and feedback experiments','The spatial envelope is (1+cos(2πx/L))/2, evaluated at each oxygen. The dipole potential includes the corresponding gradient force, so this test changes translational forcing as well as torque. No molecule ID is given a privileged role. Feedback applies E=clip(2P,0.5 V/nm) or E=clip(2[-Py,Px,0],0.5 V/nm), updated every 0.1 ps. Any collective rotation or ordering in these cases belongs to the coupled water–controller system.'),
        ('Recurrence and reduced equations','The analysis tries 2, 3, 4 and 6 coarse states without calling them saddles. It fits on baseline start 1 and evaluates starts 2 and 3. A transition matrix, PCA reduction or recurring state label does not demonstrate a heteroclinic network. Such a claim needs invariant saddle states, their stable/unstable directions, and connecting trajectories in a validated dynamical description.'),
        ('Numerical validation',f'The field-gradient and permutation checks pass. In a separate 0.2 ps deterministic audit, halving the MD step and field update gives coarse/fine error ratios {validation["MD_step_error_reduction_ratio"]:.2f} and {validation["field_update_error_reduction_ratio"]:.2f} relative to a finer reference. External work–energy residuals are recorded rather than forced to zero; the largest is {max(r["relative_work_energy_residual"] for r in validation["runs"]):.2e} of initial total energy. Twelve stochastic refinement trajectories compare ensemble summaries. They are not expected to reproduce pointwise paths.'),
        ('What the mathematical examples establish','The 12-family reference suite tests Adler entrainment, the capacitive Josephson equation, protein feedback, cusp geometry, delayed onset, attraction versus Lyapunov stability, Maxwell–Bloch reduction, SIR thresholds, identical phase oscillators, three-state heteroclinic competition, engineered nine-state groups, and spatial activity coupling. The circle example has θ′=1−cosθ: trajectories approach the same point eventually, while arbitrarily close initial points can first travel far away. Parameters and group structure in these examples are explicitly chosen. They test mechanisms and analysis behavior; they do not demonstrate those mechanisms in water.'),
        ('Scope','This is a finite exploratory battery in a small, periodic, rigid fixed-charge model. It contains no electronic polarization, chemical reactions, interfaces, electrodes or superconducting degrees of freedom. Strong applied fields are model probes. Stable molecular identities, permanent hierarchy, biological agency and quantum coherence are not inferred. A weak or absent effect is part of the result.')]
    old=[('Positional sensitivity and finite-time impulses','../response-extension/README.md'),('Drag, integration and rotating hoop','../damping-extension/README.md'),('Switching and phase models','../switching-extension/README.md'),('Linear phase portraits, stability and reversibility','../stability-extension/README.md'),('Pendulum, parametric swing, Van der Pol, Duffing, Hopf and pitchfork','../oscillator-extension/README.md'),('Static electric-field water experiment','../driven-water-extension/README.md')]
    headers=['Condition','Water |P|, mean ± SD','Control |P|, mean ± SD','Water temperature (K)','Water tracking screen']
    table_rows=[]
    for case in cases():
        a=aggregate.get(('water',case['name']));b=aggregate.get(('rotors',case['name']))
        table_rows.append([case['name'].replace('_',' '),ms(get('water',case['name'],'late_norm')),ms(get('rotors',case['name'],'late_norm')),ms(get('water',case['name'],'late_temperature_K')),(str(a['tracking_count'])+'/'+str(a['late_norm']['n']) if a else 'pending') if case['mode']=='rotate' else 'not applicable'])
    def mdtable(headers,rows):
        cell=lambda value:value.replace('|',r'\|')
        return '| '+' | '.join(map(cell,headers))+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(cell,row))+' |' for row in rows)
    commands='python exploratory-water-extension/water_battery.py\npython exploratory-water-extension/water_battery.py --refinements\npython exploratory-water-extension/model_benchmarks.py\npython exploratory-water-extension/numerical_validation.py\npython exploratory-water-extension/analyze_battery.py\npython exploratory-water-extension/build_report.py'
    title='A broad search for collective behavior in identical water'
    rotation_headers=['Frequency (cycles/ps)','Water in-phase alignment','Water coherent amplitude','Water lag (degrees)','Control in-phase alignment','Control coherent amplitude','Control lag (degrees)']
    rotation_note='All rows use a 0.5 V/nm rotating field and the same phase-analysis window within each run; entries are means over three configurations. In-phase alignment is the real part of the rotating response; its magnitude measures coherent amplitude. Lag is the angle of this averaged response, not a claim of reliable instantaneous phase when amplitude is small. The steady-field values above use the last third of each run.'
    text='# '+title+'\n\n'+('\n\n'.join(statements))+'\n\n## Following up the steady-versus-rotating contrast\n\n'+rotation_note+'\n\n'+mdtable(rotation_headers,rotation_comparison)+'\n\n## Every condition\n\n'+mdtable(headers,table_rows)
    for name,caption in figures:text+='\n\n## '+caption+'\n\n!['+caption+']('+name+'.png)'
    for heading,body in method:text+='\n\n## '+heading+'\n\n'+body
    text+='\n\n## Earlier completed tests\n\n'+'\n'.join('- ['+label+']('+url+')' for label,url in old)
    text+='\n\n## Reproduce and inspect\n\nRun from the parent experiment folder with the versions in requirements.txt. The three saved response-extension inputs and original system are required; hashes and seeds are embedded in every run archive. Resume checks reject mismatched configurations.\n\n```sh\n'+commands+'\n```\n\n[Protocol](data/protocol.json) · [All results](data/results.json) · [Run summary CSV](data/run_summary.csv) · [Model results](data/model_results.json) · [Numerical audit](data/integration_validation.json) · [Checksums](manifest_sha256.json)\n\nThe NPZ files retain 0.1 ps collective/local observables, 0.5 ps molecular orientations, five spatial snapshots and final positions/velocities. They are reduced recorded trajectories, not every integration step.\n\n## References\n\nThe supplied screenshots motivate the reference equations. [Heteroclinic networks for brain dynamics](https://www.frontiersin.org/journals/network-physiology/articles/10.3389/fnetp.2023.1276401/full) provides the network context; our coefficients and water mappings are explicitly specified experiments. [OpenMM platform documentation](https://docs.openmm.org/latest/userguide/library/04_platform_specifics.html) describes the numerical precision settings. All production runs here use OpenCL double precision.\n'
    (target/'README.md').write_text(text,encoding='utf-8')
    html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Collective water experiment battery</title><style>body{margin:0;background:#f7fafb;color:#193440;font:16px/1.65 system-ui,sans-serif}main{max-width:1180px;margin:auto;padding:42px 24px}h1{font-size:44px;line-height:1.1;max-width:920px}h2{font-size:25px;margin-top:42px}small,.eyebrow{color:#55707a}.intro{background:#e7f2f3;border-left:5px solid #087e8b;padding:20px 24px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;text-align:left;border-bottom:1px solid #ccdce0}th{background:#e7eff2}img{width:100%;height:auto;display:block;margin:20px 0}a{color:#087e8b}input{padding:10px;border:1px solid #9db7c1;border-radius:6px;width:min(500px,90%);margin:10px 0}pre{overflow:auto;background:#e8f0f2;padding:18px}.note{color:#586d76;font-size:14px}li{margin:9px 0}</style></head><body><main>'''
    html+='<p class="eyebrow">EXPLORATORY MOLECULAR STUDY · REPRODUCIBLE CONTROLS · 9 OCTOBER 2026</p><h1>'+title+'</h1><div class="intro"><strong>Try the mechanisms. Keep the outcomes.</strong><p>A systematic battery of changing drives, histories, spatial conditions and feedback, with water simulations and mathematical reference models clearly identified.</p></div>'
    if args.preview:html+='<p><strong>INCOMPLETE PREVIEW — not the final release.</strong></p>'
    html+='<h2>What the experiments found</h2><ul>'+''.join('<li>'+escape(s)+'</li>' for s in statements)+'</ul><h2>Following up the steady-versus-rotating contrast</h2><p>'+escape(rotation_note)+'</p><div class="scroll"><table><thead><tr>'+''.join('<th>'+escape(h)+'</th>' for h in rotation_headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+escape(c)+'</td>' for c in row)+'</tr>' for row in rotation_comparison)+'</tbody></table></div>'
    html+='<h2>All 26 conditions</h2><p class="note">|P| is collective orientation magnitude, not a fraction of aligned molecules. Values after ± are between-run SD. Filter the table to inspect a family.</p><input id="filter" placeholder="Filter: rotate, history, feedback, sequence…" aria-label="Filter experiment table"><div class="scroll"><table id="experiments"><thead><tr>'+''.join('<th>'+escape(h)+'</th>' for h in headers)+'</tr></thead><tbody>'
    html+=''.join('<tr>'+''.join('<td>'+escape(c)+'</td>' for c in row)+'</tr>' for row in table_rows)+'</tbody></table></div>'
    for name,caption in figures:html+='<h2>'+escape(caption)+'</h2><img alt="'+escape(caption)+'" src="data:image/png;base64,'+base64.b64encode((target/(name+'.png')).read_bytes()).decode()+'">'
    for heading,body in method:html+='<h2>'+escape(heading)+'</h2><p>'+escape(body)+'</p>'
    html+='<h2>Earlier completed tests</h2><ul>'+''.join('<li><a href="'+url+'">'+escape(label)+'</a></li>' for label,url in old)+'</ul><h2>Reproduce and inspect</h2><pre>'+escape(commands)+'</pre><p><a href="README.md">Methods and file index</a> · <a href="data/results.json">All results</a> · <a href="data/run_summary.csv">Run summary</a> · <a href="data/protocol.json">Protocol</a></p><p class="note">Reduced trajectories and checkpoints are retained. References and limitations appear in the methods. No physical water experiment was performed.</p></main><script>document.getElementById("filter").addEventListener("input",function(){let q=this.value.toLowerCase();document.querySelectorAll("#experiments tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></body></html>'
    (target/'report.html').write_text(html,encoding='utf-8')
    if not args.preview:
        files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ['manifest_sha256.json','partial_results.json']]
        save(HERE/'manifest_sha256.json',{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)})
    print(json.dumps(dict(figures=len(figures),report=str(target/'report.html'),water_rotation_screens=positives,rotor_rotation_screens=rotorpositives)),flush=True)


if __name__=='__main__':main()
