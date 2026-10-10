"""Analyze replica means without treating correlated frames as independent."""
from pathlib import Path
import json
import numpy as np
from scipy.stats import rankdata, spearmanr, t as student_t
from drive_water import save, SOURCES, TIMES, KB, C

HERE=Path(__file__).resolve().parent;D=HERE/'data'


def summarize(values):
    a=np.asarray(values,dtype=float);n=len(a)
    sd=float(np.std(a,ddof=1)) if n>1 else 0.
    half=float(student_t.ppf(.975,n-1)*sd/np.sqrt(n)) if n>1 else None
    return {'replica_values':a.tolist(),'n_replicas':n,'mean':float(a.mean()),'SD_between_replicas':sd,
            'approximate_95pct_t_interval':[float(a.mean()-half),float(a.mean()+half)] if half is not None else None}


def run():
    protocol=json.loads((D/'protocol.json').read_text());reps=protocol['replicas']
    assert len(reps)==3,'The delivered design requires all three initial states.'
    rows=[];loaded={}
    for kind,fields in [('water',[0.,.1,.5,-.5]),('rotors',[0.,.5])]:
        for rep in reps:
            for field in fields:
                name=f'{kind}_r{rep}_E{field:g}'
                data=np.load(D/(name+'.npz'));loaded[kind,rep,field]=data
                assert all(np.isfinite(data[key]).all() for key in data.files)
                t=data['time_ps'];on=(t>=6-1e-8)&(t<=10+1e-8);off=(t>=16-1e-8)&(t<=20+1e-8)
                p=np.asarray(data['dipole_directions'],float).mean(axis=1)
                coord=data['coordination_scores'];hb=data['hbond_counts']
                row={'kind':kind,'replica':rep,'field_V_nm':field,'on_alignment_z':float(p[on,2].mean()),
                     'recovery_alignment_z':float(p[off,2].mean()),
                     'on_temperature_K':float(data['temperature_K'][on].mean()),
                     'recovery_temperature_K':float(data['temperature_K'][off].mean()),
                     'on_connected_neighbor_orientation':float(data['connected_neighbor_orientation'][on].mean()),
                     'on_hbond_neighbors':float(hb[on].mean()),
                     'on_coordination_cv':float(np.mean(coord[on].std(axis=1)/coord[on].mean(axis=1)))}
                if kind=='water':
                    ranks=rankdata(coord[t<=10+1e-8],axis=1);ranks-=ranks.mean(axis=1,keepdims=True)
                    row['coordination_rank_correlation_by_lag_ps']={}
                    for lag in [1,10,50]:
                        a,b=ranks[:-lag],ranks[lag:]
                        value=np.mean(np.sum(a*b,axis=1)/np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1)))
                        row['coordination_rank_correlation_by_lag_ps'][str(lag*.1)]=float(value)
                    if field>=0:
                        score=data['impulse_response_norms_ps'].sum(axis=2)
                        row['selected_source_mean_scores_ps']=score.mean(axis=0).tolist()
                        row['selected_source_scores_200fs_ps']=score[:,-1].tolist()
                    if 'force_scores' in data:
                        force=data['force_scores']
                        row.update({'snapshot_force_leader':int(np.argmax(force)),
                                    'snapshot_force_mean':float(force.mean()),
                                    'snapshot_force_cv':float(force.std()/force.mean()),
                                    'snapshot_force_max_min':float(force.max()/force.min()),
                                    'force_jacobian_reciprocity_error':float(data['force_jacobian_reciprocity_error'])})
                rows.append(row)
    def values(kind,field,key):return [row[key] for row in rows if row['kind']==kind and row['field_V_nm']==field]
    aggregates=[]
    for kind,fields in [('water',[0.,.1,.5,-.5]),('rotors',[0.,.5])]:
        for field in fields:
            aggregates.append({'kind':kind,'field_V_nm':field,
              **{key:summarize(values(kind,field,key)) for key in ['on_alignment_z','recovery_alignment_z','on_temperature_K','on_connected_neighbor_orientation','on_hbond_neighbors']}})
    p0=np.array(values('water',0.,'on_alignment_z'));p1=np.array(values('water',.1,'on_alignment_z'));p5=np.array(values('water',.5,'on_alignment_z'))
    differences={
        'alignment_gain_0p1_vs_zero':summarize(p1-p0),
        'alignment_gain_0p5_vs_zero':summarize(p5-p0),
        'departure_from_fivefold_linear_scaling':summarize((p5-p0)-5*(p1-p0)),
        'recovery_residual_0p5_vs_zero':summarize(np.array(values('water',.5,'recovery_alignment_z'))-values('water',0.,'recovery_alignment_z')),
        'opposite_fields_alignment_sum':summarize(p5+np.array(values('water',-.5,'on_alignment_z'))),
        'interacting_minus_noninteracting_alignment_0p5':summarize(p5-np.array(values('rotors',.5,'on_alignment_z'))),
        'hbond_neighbor_change_0p5_vs_zero':summarize(np.array(values('water',.5,'on_hbond_neighbors'))-values('water',0.,'on_hbond_neighbors'))}
    response=[]
    for field in [.1,.5]:
        ratios=[]
        for rep in reps:
            a=loaded['water',rep,0.]['impulse_response_norms_ps'].sum(axis=2)
            b=loaded['water',rep,field]['impulse_response_norms_ps'].sum(axis=2)
            ratios.append(float(b[:,-1].mean()/a[:,-1].mean()))
        response.append({'field_V_nm':field,'ratio_of_six_source_mean_200fs_response':summarize(ratios)})
    controlled=loaded['water',0,.5]
    base=controlled['impulse_response_norms_ps'].sum(axis=2)
    def relative(key):return (np.linalg.norm(controlled[key]-base,axis=0)/np.linalg.norm(base,axis=0)).tolist()
    controls=json.loads((D/'controls.json').read_text())
    controls.update({'half_step_relative_L2_by_time':relative('half_step_scores_ps'),
                    'half_kick_relative_L2_by_time':relative('half_kick_scores_ps'),
                    'permutation_relative_L2_by_time':relative('permuted_impulse_scores_ps'),
                    'noninteracting_max_impulse_score_ps':float(controlled['free_impulse_scores_ps'].max())})
    force=controlled['force_scores'];perm=controlled['permutation'];pforce=controlled['permuted_force_scores']
    controls['force_score_permutation_relative_L2']=float(np.linalg.norm(pforce-force[perm])/np.linalg.norm(force))
    controls['original_force_leader']=int(force.argmax())
    controls['permuted_force_leader']=int(pforce.argmax())
    controls['expected_permuted_force_leader']=int(np.flatnonzero(perm==force.argmax())[0])
    controls['acceptance']={'half_step_200fs_relative_L2_below':.02,'half_kick_max_relative_L2_below':.02,
                            'permutation_relative_L2_below':1e-5,'noninteracting_score_below_ps':1e-7}
    assert controls['half_step_relative_L2_by_time'][-1]<.02
    assert max(controls['half_kick_relative_L2_by_time'])<.02
    assert max(controls['permutation_relative_L2_by_time'])<1e-5
    assert controls['force_score_permutation_relative_L2']<1e-5
    assert controls['noninteracting_max_impulse_score_ps']<1e-7
    assert controls['permuted_force_leader']==controls['expected_permuted_force_leader']
    controls['passed']=True
    # Independent single-dipole equilibrium prediction; these runs remain finite-duration checks.
    x=loaded['rotors',0,.5]['snapshot_positions_nm'][0]
    side=protocol['box_side_nm'];oh=x[:,1:]-x[:,0:1];oh-=side*np.rint(oh/side)
    dipole=float(np.linalg.norm(.417*oh.sum(axis=1),axis=1).mean())
    alpha=C*.5*dipole/(KB*300)
    predicted=float(1/np.tanh(alpha)-1/alpha)
    # Preserve each frame's orientation distribution and geometry separately;
    # shuffle only their association. This is a statistical null, not dynamics.
    rng=np.random.default_rng(20261400);spatial=[]
    side=protocol['box_side_nm']
    for kind in ['water','rotors']:
        for rep in reps:
            for field in [0.,.5]:
                data=loaded[kind,rep,field]
                for snapshot,time_ps in enumerate([0.,10.,20.]):
                    xyz=data['snapshot_positions_nm'][snapshot]
                    oh=xyz[:,1:]-xyz[:,0:1];oh-=side*np.rint(oh/side)
                    u=oh.sum(axis=1);u/=np.linalg.norm(u,axis=1)[:,None]
                    oo=xyz[:,0,None,:]-xyz[None,:,0,:];oo-=side*np.rint(oo/side)
                    dist=np.linalg.norm(oo,axis=2)
                    ii,jj=np.where(np.triu((dist<.35)&(dist>1e-8),1))
                    dots=u@u.T;observed=float(dots[ii,jj].mean())
                    null=[]
                    for _ in range(200):
                        permutation=rng.permutation(len(u));null.append(float(dots[permutation[ii],permutation[jj]].mean()))
                    exact_null=float((len(u)*np.dot(u.mean(axis=0),u.mean(axis=0))-1)/(len(u)-1))
                    spatial.append({'kind':kind,'replica':rep,'field_V_nm':field,'time_ps':time_ps,
                      'neighbor_pairs':len(ii),'observed_neighbor_dot':observed,
                      'exact_shuffled_expectation':exact_null,'excess_over_shuffled_expectation':observed-exact_null,
                      'shuffle_2p5_97p5_percentiles':np.percentile(null,[2.5,97.5]).tolist(),
                      'shuffle_mean':float(np.mean(null)),'shuffle_SD':float(np.std(null,ddof=1)),
                      'observed_above_shuffle_97p5':observed>np.percentile(null,97.5)})
    relaxation=[]
    for kind in ['water','rotors']:
        for rep in reps:
            data=loaded[kind,rep,.5];time=data['time_ps'];mask=time>=10-1e-8
            elapsed=time[mask]-10;alignment=data['dipole_directions'][mask,:,2].mean(axis=1)
            level=float(alignment[0]/2);hits=np.flatnonzero(alignment<=level)
            k=int(hits[0]) if len(hits) else None
            bracket=[float(elapsed[k-1]),float(elapsed[k])] if k is not None and k>0 else None
            relaxation.append({'kind':kind,'replica':rep,'initial_field_off_alignment':float(alignment[0]),
               'half_alignment_threshold':level,'first_half_decay_bracket_ps':bracket,
               'interpretation':'First observed drop below half of alignment at field removal, bracketed by saved frames. No baseline subtraction, fitted exponential, or equilibrium claim.'})
    result={'rows':rows,'aggregates':aggregates,'contrasts':differences,'conditional_response':response,
            'field_off_half_decay':relaxation,
            'spatial_orientation_test':{'cutoff_nm':.35,'shuffle_count':200,'seed':20261400,'rows':spatial,
             'null':'Permute the observed dipole directions across fixed oxygen positions, preserving the complete global orientation distribution. These shuffled configurations are not simulated physical states.',
             'meaning':'Neighbor alignment above this null is a spatial association beyond common global alignment; it alone does not establish intent, leadership, or long-lived domains.'},
            'noninteracting_equilibrium_prediction':{'dipole_e_nm':dipole,'alpha':alpha,'Langevin_alignment':predicted,
              'formula':'L(alpha)=coth(alpha)-1/alpha; alpha=96.48533212331002*E*p/(R*T)'},
            'controls':controls,'passed':True,
            'interpretation_rules':['Alignment from a common external field does not by itself establish cooperation.',
              'Six selected impulse sources cannot determine the global finite-time leader.',
              'A nonzero short-time recovery residual is not evidence of permanent memory or hysteresis.',
              'Comparing three replica means is exploratory; reported t intervals assume approximately normal replica means and do not replace larger ensembles.']}
    save(D/'results.json',result)
    print(json.dumps({'contrasts':differences,'response':response,'controls':controls},indent=2),flush=True)


if __name__=='__main__':run()
