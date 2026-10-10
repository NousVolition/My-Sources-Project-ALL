"""Render verified driven-water measurements and their interpretation."""
from pathlib import Path
from html import escape
import base64,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;D=HERE/'data'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.titleweight':'bold','figure.facecolor':'white',
                     'axes.prop_cycle':plt.cycler(color=['#78858d','#3e87c2','#087e83','#df7636'])})


def main():
    r=json.loads((D/'results.json').read_text());p=json.loads((D/'protocol.json').read_text())
    assert r['passed'] and r['controls']['passed']
    def group(kind,field):return next(a for a in r['aggregates'] if a['kind']==kind and a['field_V_nm']==field)
    def load(kind,rep,field):return np.load(D/f'{kind}_r{rep}_E{field:g}.npz')
    def save(fig,name):
        fig.tight_layout(pad=2.);fig.savefig(HERE/(name+'.png'),dpi=150,bbox_inches='tight')
        fig.savefig(HERE/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(12,9))
    for field in [0.,.1,.5,-.5]:
        trials=[load('water',rep,field) for rep in range(3)];t=trials[0]['time_ps']
        trace=np.array([a['dipole_directions'][:,:,2].mean(axis=1) for a in trials])
        mean=trace.mean(axis=0);sd=trace.std(axis=0,ddof=1)
        line=ax[0,0].plot(t,mean,label=f'{field:g} V/nm')[0]
        ax[0,0].fill_between(t,mean-sd,mean+sd,color=line.get_color(),alpha=.1)
    ax[0,0].axvline(10,color='black',ls=':',lw=1)
    ax[0,0].set(title='Switching the field changes water orientation',xlabel='Time (ps); field removed at 10 ps',ylabel='Mean dipole direction along z')
    ax[0,0].legend(fontsize=9);ax[0,0].text(.03,.05,'shading: SD across 3 runs',transform=ax[0,0].transAxes,fontsize=9)
    for kind,label in [('water','interacting model water'),('rotors','noninteracting rigid dipoles')]:
        trials=[load(kind,rep,.5) for rep in range(3)]
        trace=np.array([a['dipole_directions'][:,:,2].mean(axis=1) for a in trials])
        ax[0,1].plot(t,trace.mean(axis=0),label=label)
    prediction=r['noninteracting_equilibrium_prediction']['Langevin_alignment']
    ax[0,1].hlines(prediction,0,10,color='black',ls='--',lw=1,label='independent-dipole equilibrium prediction')
    ax[0,1].axvline(10,color='black',ls=':',lw=1)
    ax[0,1].set(title='Shared forcing can align noninteracting dipoles',xlabel='Time (ps)',ylabel='Mean dipole direction along z');ax[0,1].legend(fontsize=8)
    fields=np.array([0.,.1,.5]);means=np.array([group('water',f)['on_alignment_z']['mean'] for f in fields])
    sd=np.array([group('water',f)['on_alignment_z']['SD_between_replicas'] for f in fields])
    ax[1,0].errorbar(fields,means,yerr=sd,fmt='o-',label='water: mean ± replica SD')
    ax[1,0].plot(fields,means[0]+(means[1]-means[0])*fields/.1,'--',label='linear extrapolation from 0 → 0.1')
    ax[1,0].set(title='Test proportionality instead of assuming it',xlabel='Applied field (V/nm)',ylabel='Mean alignment during 6–10 ps');ax[1,0].legend(fontsize=8)
    for kind,offset in [('water',-.03),('rotors',.03)]:
        values=[group(kind,f)['on_alignment_z'] for f in [0.,.5]]
        ax[1,1].errorbar(np.array([0.,.5])+offset,[v['mean'] for v in values],yerr=[v['SD_between_replicas'] for v in values],fmt='o',label=kind)
    ax[1,1].set(title='Intermolecular interactions change the response',xlabel='Applied field (V/nm); slight plotting offset',ylabel='Mean alignment during 6–10 ps');ax[1,1].legend()
    save(fig,'collective_response')

    fig,ax=plt.subplots(2,2,figsize=(12,9))
    for field in [0.,.5]:
        scores=np.concatenate([load('water',rep,field)['force_scores'] for rep in range(3)])
        ax[0,0].plot(np.sort(scores),np.arange(1,len(scores)+1)/len(scores),label=f'{field:g} V/nm')
    ax[0,0].set(title='Snapshot force influence remains unequal',xlabel='Sum of off-source force-sensitivity norms (kJ mol⁻¹ nm⁻²)',ylabel='Fraction of pooled molecule scores');ax[0,0].legend()
    for field in [.1,.5]:
        values=next(a['ratio_of_six_source_mean_200fs_response'] for a in r['conditional_response'] if a['field_V_nm']==field)
        ax[0,1].scatter(np.full(3,field),values['replica_values'],s=35,label=f'{field:g} V/nm')
        ax[0,1].plot([field-.025,field+.025],[values['mean']]*2,color='black',lw=2)
    ax[0,1].axhline(1,color='gray',ls='--')
    ax[0,1].set(title='Six prespecified molecules: finite-time response',xlabel='Field maintained during each probe (V/nm)',ylabel='Mean 200-fs response / matched zero-field mean')
    for field in [0.,.5]:
        rows=[a for a in r['rows'] if a['kind']=='water' and a['field_V_nm']==field]
        lags=[.1,1.,5.]
        vals=np.array([[a['coordination_rank_correlation_by_lag_ps'][str(lag)] for lag in lags] for a in rows])
        ax[1,0].errorbar(lags,vals.mean(axis=0),yerr=vals.std(axis=0,ddof=1),fmt='o-',label=f'{field:g} V/nm')
    ax[1,0].set_xscale('log');ax[1,0].set(title='Local-neighborhood ranks change over time',xlabel='Lag (ps)',ylabel='Rank correlation, mean ± replica SD');ax[1,0].legend()
    for field in [0.,.1,.5,-.5]:
        trials=[load('water',rep,field) for rep in range(3)]
        ax[1,1].plot(t,np.mean([a['temperature_K'] for a in trials],axis=0),lw=.7,label=f'{field:g} V/nm')
    ax[1,1].axhline(300,color='black',ls='--',lw=1)
    ax[1,1].set(title='Thermal bath maintained at 300 K',xlabel='Time (ps)',ylabel='Instantaneous kinetic temperature (K)');ax[1,1].legend(fontsize=8)
    save(fig,'influence_controls')

    fig,ax=plt.subplots(1,2,figsize=(12,4.5))
    spatial=r['spatial_orientation_test']['rows']
    for column,field in enumerate([0.,.5]):
        for kind,offset,color in [('water',-.15,'#087e83'),('rotors',.15,'#df7636')]:
            records=[a for a in spatial if a['kind']==kind and a['field_V_nm']==field and a['time_ps']==10.]
            obs=np.array([a['observed_neighbor_dot'] for a in records]);null=np.array([a['exact_shuffled_expectation'] for a in records])
            lower=np.array([a['shuffle_2p5_97p5_percentiles'][0] for a in records]);upper=np.array([a['shuffle_2p5_97p5_percentiles'][1] for a in records])
            xx=np.arange(3)+offset
            ax[column].vlines(xx,lower-null,upper-null,color=color,lw=5,alpha=.25,label=kind+' shuffle 95% range')
            ax[column].scatter(xx,obs-null,color=color,s=45,label=kind+' observed excess')
        ax[column].axhline(0,color='black',lw=.7,ls='--')
        ax[column].set(title=f'Local orientation pattern at E = {field:g} V/nm',xlabel='Replica',ylabel='Neighbor alignment minus shuffled expectation',xticks=[0,1,2])
        ax[column].legend(fontsize=8)
    save(fig,'spatial_patterns')

    fig,ax=plt.subplots(1,3,figsize=(13,4.7))
    trial=load('water',0,.5);side=p['box_side_nm']
    for i,title in enumerate(['Initial state','After 10 ps at +0.5 V/nm','After 10 ps with field removed']):
        xyz=trial['snapshot_positions_nm'][i];oxygen=xyz[:,0]%side
        local=xyz[:,1:]-xyz[:,0:1];local-=side*np.rint(local/side)
        direction=local.sum(axis=1);direction/=np.linalg.norm(direction,axis=1)[:,None]
        scatter=ax[i].scatter(oxygen[:,0],oxygen[:,2],c=direction[:,2],cmap='coolwarm',vmin=-1,vmax=1,s=18,alpha=.85)
        ax[i].quiver(oxygen[:,0],oxygen[:,2],direction[:,0],direction[:,2],angles='xy',scale_units='xy',scale=14,width=.003,alpha=.6)
        ax[i].set(title=title,xlabel='x position (nm)',ylabel='z position (nm)',xlim=(0,side),ylim=(0,side));ax[i].set_aspect('equal')
    fig.tight_layout(rect=[0,0,1,.9]);fig.suptitle('Actual saved molecular snapshots: arrows show dipole directions',y=.99,fontsize=15)
    # Separate legend text avoids squeezing the axes with an overlapping colorbar.
    fig.text(.5,.04,'Projection onto x–z plane; red = dipole points toward +z, blue = toward −z. Arrow length is a display scale.',ha='center',fontsize=9)
    fig.subplots_adjust(bottom=.2,top=.82)
    fig.savefig(HERE/'water_snapshots.png',dpi=150,bbox_inches='tight');fig.savefig(HERE/'water_snapshots.svg',bbox_inches='tight');plt.close(fig)

    a0=group('water',0.);a1=group('water',.1);a5=group('water',.5);aN=group('rotors',.5)
    contrasts=r['contrasts'];checks=r['controls']
    fmt=lambda value:f'{value["mean"]:.4f} ± {value["SD_between_replicas"]:.4f}'
    table_rows=[[f'{f:g}',fmt(group('water',f)['on_alignment_z']),fmt(group('water',f)['recovery_alignment_z']),fmt(group('water',f)['on_temperature_K'])] for f in [0.,.1,.5,-.5]]
    difference=contrasts['departure_from_fivefold_linear_scaling'];ci=difference['approximate_95pct_t_interval']
    scale_conclusion=('The replica-mean contrast and its exploratory interval are below zero, supporting a sublinear response over these tested fields.' if ci[1]<0 else 'The exploratory interval does not resolve a sublinear response; more replicas are needed.' )
    strong_spatial=[a for a in spatial if a['kind']=='water' and a['field_V_nm']==.5 and a['time_ps']==10.]
    spatial_excess=np.array([a['excess_over_shuffled_expectation'] for a in strong_spatial])
    def decay_brackets(kind):
        return ', '.join('not reached' if a['first_half_decay_bracket_ps'] is None else f'{a["first_half_decay_bracket_ps"][0]:.1f}–{a["first_half_decay_bracket_ps"][1]:.1f} ps' for a in r['field_off_half_decay'] if a['kind']==kind)
    paragraphs=[
      ('What changed','This experiment changes the conditions experienced by actual simulated water molecules. All 216 molecules retain identical masses, rigid geometry, charges and pair interactions. A uniform electric field couples to those charges for 10 ps and is removed for the next 10 ps. We compare 0, +0.1, +0.5 and −0.5 V/nm, from three previously saved water states each, with a 300 K Langevin bath. All trajectories are retained.'),
      ('What alignment means','Let uᵢ be the unit vector along each molecular dipole and Pz=mean(uᵢ,z). Pz ranges from −1 to +1. It describes average orientation, not a fraction of molecules, a leader score, or intent. Window means during 6–10 ps are '+fmt(a0['on_alignment_z'])+' at zero field and '+fmt(a5['on_alignment_z'])+' at +0.5 V/nm. Values after ± are sample standard deviations across three replica means.'),
      ('A direct nonlinearity test','If the mean field-induced alignment were proportional to the applied field, ΔP(0.5) would equal 5ΔP(0.1), where ΔP(E)=P(E)−P(0). The measured contrast ΔP(0.5)−5ΔP(0.1) is '+fmt(difference)+f'; its exploratory 95% t interval is [{ci[0]:.4f}, {ci[1]:.4f}]. '+scale_conclusion+' This tests the observed response curve. It does not make every individual molecular interaction proportional at weak field.'),
      ('Shared forcing versus interaction','Noninteracting rigid dipoles also align: at +0.5 V/nm their window mean is '+fmt(aN['on_alignment_z'])+'. The independent-dipole equilibrium prediction is '+f'{prediction:.4f}'+'. The interacting-minus-noninteracting difference is '+fmt(contrasts['interacting_minus_noninteracting_alignment_0p5'])+'. This distinguishes the shared external stimulus from the effect of intermolecular interactions. The noninteracting control contains overlapping ghost particles and is not liquid water. Alignment alone is insufficient evidence of cooperation.'),
      ('Do neighbors form an additional orientation pattern?','At each saved initial, field-on-end and recovery-end snapshot, we hold oxygen positions fixed and randomly reassign the existing dipole directions 200 times. This preserves the complete overall orientation distribution. For O–O pairs closer than 0.35 nm, we compare the observed mean dipole dot product with that spatially shuffled null. At the +0.5 V/nm endpoint, the three interacting-water excesses are '+', '.join(f'{value:.4f}' for value in spatial_excess)+'. '+str(sum(a['observed_above_shuffle_97p5'] for a in strong_spatial))+' of three exceed their own 97.5th-percentile shuffle value. The plot also shows zero-field and noninteracting controls. This tests a spatial association beyond shared global alignment; the shuffled states are a statistical reference, not physical water trajectories. It does not demonstrate intent or a long-lived collective domain.'),
      ('Does the change persist?','The field is turned off at 10 ps. During 16–20 ps, the +0.5-minus-zero-field alignment difference is '+fmt(contrasts['recovery_residual_0p5_vs_zero'])+'. This is a measured recovery-window residual. Twenty picoseconds cannot establish permanent memory, metastability, or hysteresis. A memory claim would require much longer field-free runs, repeated histories, and comparison of distinct stable states under the same final conditions.'),
      ('Influence without privileged identities','At the 10 ps endpoints for E=0 and E=0.5, we calculate the full translational force Jacobian for all 216 molecules, as in the original study. Snapshot influence is the sum of off-source 3×3 block norms. It remains a reciprocal, geometry-dependent sensitivity, not a directional hierarchy. Permuting complete molecular states moves the leading label from '+str(checks['original_force_leader'])+' to '+str(checks['permuted_force_leader'])+', exactly the expected relabeling; relative score error is '+f'{checks["force_score_permutation_relative_L2"]:.2e}'+'. Leader changes between separately evolved trajectories also occur from ordinary molecular motion, so such changes alone cannot be attributed to the field.'),
      ('Finite-time intervention, with a limited scope','We probe source IDs 0,43,86,129,172,215, chosen before seeing outcomes, using paired ±0.01 nm/ps kicks along all three axes. Each source kick is balanced by compensating kicks to the other molecules. Target center-of-mass response is corrected for that direct compensation; score is summed off-source response-block norms in ps. The field remains on during these deterministic 0.2 ps probes, with no thermostat. The source positions are taken from the 10 ps endpoint; fresh constrained Maxwell velocities at 300 K define each conditional probe. At +0.5 V/nm the mean ratio of six-source 200-fs response to the matching zero-field mean is '+fmt(r['conditional_response'][1]['ratio_of_six_source_mean_200fs_response'])+'. Six sources do not identify the global finite-time leader, and this short-horizon ratio is not an energy-transfer fraction.'),
      ('Numerical controls','The field force agrees with qE to '+f'{checks["force_max_error_kJ_mol_nm"]:.2e}'+ ' kJ mol⁻¹ nm⁻¹. A periodic-image shift changes its energy by '+f'{checks["periodic_image_energy_difference_kJ_mol"]:.2e}'+ ' kJ/mol. On the selected +0.5 V/nm replica, halving the probe step changes 200-fs scores by '+f'{100*checks["half_step_relative_L2_by_time"][-1]:.4f}%'+ ' in relative L2 norm; the largest half-kick difference over measured times is '+f'{100*max(checks["half_kick_relative_L2_by_time"]):.4f}%'+'. The no-interaction impulse score is at most '+f'{checks["noninteracting_max_impulse_score_ps"]:.2e}'+ ' ps. Permutation, finite-array and stated acceptance checks pass. These deterministic checks do not substitute for full stochastic-ensemble time-step convergence.'),
      ('Methods and uncertainty','The existing OpenMM TIP3P model is rigid and nonpolarizable, in a cubic periodic box of side '+f'{p["box_side_nm"]:.6f}'+ ' nm at fixed density 0.997 g/cm³. Sampling uses 1 fs steps, friction 1 ps⁻¹, constraint tolerance 10⁻⁸ and 0.1 ps saved observations. Each condition starts from the same saved state within a replica but uses a separately specified thermostat seed. Frames and molecules are correlated; each replica window mean is one observation. Three initial states provide an exploratory ensemble. The t intervals assume approximately normal replica means; no independent-frame p-values are reported. Hydrogen-bond counts use O–O <0.35 nm and donor O–H within 30° of donor-to-acceptor O–O. Connected neighbor orientation subtracts squared global mean orientation from the pair-average dot product for O–O <0.35 nm; it is descriptive, not a causal cooperation score.'),
      ('What this says about water','Within this molecular model, changing the environment changes the collective behavior of otherwise identical water molecules. This is a water-specific test rather than a borrowed oscillator example. The conclusion remains bounded by the model: fixed charges cannot respond electronically, rigid molecules cannot vibrate internally or break bonds, and the box has no electrodes, interfaces, ions or chemical reactions. The applied fields are strong model probes, not a claim of experimentally validated behavior at those field strengths. These results do not assign molecules needs, decisions, or permanent roles.'),
      ('Next validation','Increase the replica count, box size and observation duration; repeat with a second water model and a polarizable model; vary bath friction and integration step; map weaker fields and repeated on/off histories; and probe every molecule if a complete finite-time ranking is needed. Persistent-state switching would require showing two reproducible states under the same final conditions, beyond the direct field-imposed alignment.')]
    paragraphs.append(('Relation to the supplied heteroclinic-network diagram','The supplied caption describes nine saddle states connected by heteroclinic paths. In a water interpretation, each node would represent a collective configuration of the entire system, rather than one molecule. Our driven trajectories and spatial orientation correlations do not identify saddle invariant states or their connecting manifolds. Repeated switching under one fixed set of conditions, well-defined collective coordinates, stability analysis around candidate states, and transition-path evidence would be needed to test that hypothesis. An on/off response created by changing the external field is not evidence of a heteroclinic network.'))
    paragraphs.append(('Relation to the supplied oscillator-network equation','The supplied equation is dθᵢ/dt = ω + ε Σⱼ wᵢⱼ H(θᵢ−θⱼ). Identical ω makes the uncoupled oscillators identical; collective patterns depend on the interaction weights wᵢⱼ, coupling strength ε, function H and initial phases. The retrieved source article supplies H(θ)=sin(θ−α)−r sin(2θ). This equation provides a possible reduced description, but using it for water requires defining a measurable molecular or collective phase and deriving or fitting its interactions from molecular dynamics. We have not assigned an arbitrary intrinsic clock to every water molecule or claimed that this equation alone produces the pictured heteroclinic network. A future phase model should predict held-out molecular trajectories and retain an explicit comparison with zero coupling.'))
    paragraphs.append(('Chunking and the new reaction–diffusion equation','The newest excerpts motivate a stronger target: slow switching between groups with faster changes inside each group. A reduced model can summarize many molecular coordinates by fewer collective variables, but a useful reduction must preserve predictive behavior. The supplied candidate is ∂t Aᵢ=δ∇²Aᵢ+σAᵢ−γAᵢ²−Σⱼ≠ᵢρᵢⱼAᵢAⱼ+η|ξᵢ(t)|. Its activity variables, competition matrix, diffusion, and noise are not molecular forces. The cited Figure 5 assigns one unit a pacemaker role, so it cannot by itself demonstrate leadership emerging among equivalent units. Our present water drive is uniform and assigns no molecular leader. The accompanying NEXT_TEST.md sets out how to test nested collective patterns and any proposed reduced model without imposing the desired result.'))
    paragraphs.append(('Field removal as a relaxation test','Motivated by the supplied quench discussion, we measure the first drop below half of Pz at field removal for each +0.5 V/nm trajectory. The first crossing falls in '+decay_brackets('water')+' for the three interacting-water runs, and '+decay_brackets('rotors')+' for the noninteracting controls. Brackets come directly from 0.1 ps sampling; we do not interpolate beyond that resolution or assume exponential decay. This shows the observed relaxation timescale of an orientational response. It does not establish a bifurcation quench, underdamped oscillations, learning, or permanent storage. The latest learning excerpt motivates a separate training-then-free-recall protocol, specified in NEXT_TEST.md; that proposed learning test has not been run.'))
    def mdtable(headers,rows):return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(str(x).replace('|','\\|') for x in row)+' |' for row in rows)
    def table(headers,rows):return '<div class="scroll"><table><tr>'+''.join('<th>'+escape(x)+'</th>' for x in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</table></div>'
    headers=['Field (V/nm)','On-window Pz (mean ± SD)','Recovery-window Pz','On-window temperature (K)']
    commands='''python -m pip install -r driven-water-extension/requirements.txt
python driven-water-extension/drive_water.py --platform OpenCL
python driven-water-extension/analyze_water.py
python driven-water-extension/build_report.py'''
    sources=[('Heteroclinic networks for brain dynamics — source of the supplied network figures','https://www.frontiersin.org/journals/network-physiology/articles/10.3389/fnetp.2023.1276401/full'),
             ('OpenMM CustomCompoundBondForce','https://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomCompoundBondForce.html'),
             ('OpenMM LangevinMiddleIntegrator','https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LangevinMiddleIntegrator.html'),
             ('OpenMM water-model setup','https://docs.openmm.org/latest/userguide/application/03_model_building_editing.html')]
    readme='# Making identical model water respond differently\n\nA reproducible electric-field pulse test with interaction, thermal, numerical, and identity controls.\n\n'+mdtable(headers,table_rows)
    placements={1:'collective_response',4:'spatial_patterns',6:'influence_controls',10:'water_snapshots'}
    for i,(title,body) in enumerate(paragraphs):
        readme+='\n\n## '+title+'\n\n'+body
        if i in placements:readme+='\n\n!['+title+']('+placements[i]+'.png)'
    readme+='\n\n## Reproduce\n\nPython 3.12; run from the parent experiment folder. Requires the earlier saved response-extension/data files and data/md_protocol.json. Input hashes are recorded in the protocol.\n\n```sh\n'+commands+'\n```\n\nUse `--platform Reference` for a slower portable calculation. `--out PATH` saves trajectories elsewhere; analysis/builders read this extension’s data folder. The simulation skips existing run archives for resumption; use a fresh output directory for new parameter choices.\n\n[Report](report.html) · [Protocol](data/protocol.json) · [All numerical results and controls](data/results.json) · [Checksums](manifest_sha256.json)\n\n## Implementation references\n\n'+ '\n'.join('- ['+name+']('+url+')' for name,url in sources)+'\n\nThe custom potential uses documented coordinate expressions and global parameters; the thermal bath uses OpenMM’s documented Langevin integrator. Field and dipole formulas are derived explicitly in the protocol and code.\n\n[Original water study](../README.md) · [Earlier oscillator tests](../oscillator-extension/README.md)\n'
    (HERE/'README.md').write_text(readme,encoding='utf-8')
    html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Changing model water behavior</title><style>body{font:17px/1.65 system-ui,sans-serif;max-width:1140px;margin:40px auto;padding:0 24px;background:#fafcfd;color:#17313b}h1{font-size:42px;line-height:1.15}h2{margin-top:38px;font-size:25px}.box{background:#e7f3f3;border-left:4px solid #087e83;padding:20px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:10px;border-bottom:1px solid #d3e1e4}th{background:#e8f0f3}img{width:100%;height:auto;margin:20px 0}pre{padding:18px;background:#eaf1f4;overflow:auto}a{color:#087e83}</style></head><body><p>REPRODUCIBLE MOLECULAR SIMULATION · 9 OCTOBER 2026</p><h1>Same molecules.<br>Different collective behavior.</h1><p class="box">Change the environment, retain identical molecular properties, and test the response against undriven and noninteracting controls. These results come from the water model itself.</p>'''
    html+=table(headers,table_rows)
    for i,(title,body) in enumerate(paragraphs):
        html+='<h2>'+escape(title)+'</h2><p>'+escape(body)+'</p>'
        if i in placements:
            name=placements[i];encoded=base64.b64encode((HERE/(name+'.png')).read_bytes()).decode()
            html+='<img alt="'+escape(title)+'" src="data:image/png;base64,'+encoded+'">'
    html+='<h2>Reproduce</h2><pre>'+escape(commands)+'</pre><p>Inputs, complete numerical results, seeds and acceptance checks are documented in <a href="README.md">the methods and file index</a>.</p><h2>Implementation references</h2><p>'+ ' · '.join('<a href="'+url+'">'+escape(name)+'</a>' for name,url in sources)+'</p></body></html>'
    (HERE/'report.html').write_text(html,encoding='utf-8')
    files=[f for f in HERE.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name!='manifest_sha256.json']
    manifest={f.relative_to(HERE).as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(files)}
    save=lambda path,value:path.write_text(json.dumps(value,indent=2),encoding='utf-8')
    save(HERE/'manifest_sha256.json',manifest)
    print(json.dumps({'files':len(files)+1,'all_checks_passed':True,'alignment_zero':a0['on_alignment_z'],'alignment_0p5':a5['on_alignment_z'],'nonlinearity':difference},indent=2))


if __name__=='__main__':main()
