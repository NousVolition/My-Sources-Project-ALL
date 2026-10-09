"""Recompute tables, convergence checks, and figures from saved response arrays."""
from pathlib import Path
import json, hashlib
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent
D=R/'data'
times=np.array([.005,.02,.1,.2])
colors=['#087e83','#cf6829','#4058a5']


def native(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)


def relative_by_time(base,other):
    return {'relative_L2':np.linalg.norm(other-base,axis=0)/np.linalg.norm(base,axis=0),
            'maximum_relative':np.max(abs(other-base)/np.maximum(base,1e-30),axis=0)}


def main():
    diagnostics=[];rank_changes=[]
    fig,ax=plt.subplots(2,2,figsize=(11.7,8),layout='constrained')
    fig.suptitle('Does the snapshot leader still lead after an impulse?',fontsize=18,fontweight='bold')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    for rep in range(3):
        data=np.load(D/f'response_replica_{rep}.npz')
        scores=data['scores_ps'];force=data['force_scores'];pairs=data['pair_response_norms_ps']
        assert scores.shape==(216,4)
        assert np.allclose(scores,pairs.sum(axis=2),rtol=1e-14,atol=1e-14)
        assert np.isfinite(pairs).all() and (pairs>=0).all()
        assert np.all(pairs[np.arange(216),:,np.arange(216)]==0)
        correlations=[spearmanr(force,scores[:,k]).statistic for k in range(4)]
        ranks=[int(np.flatnonzero(np.argsort(-scores[:,k])==force.argmax())[0]+1) for k in range(4)]
        asymmetry=[np.linalg.norm(pairs[:,k,:]-pairs[:,k,:].T)/np.linalg.norm(pairs[:,k,:]) for k in range(4)]
        incoming_outgoing=[spearmanr(pairs[:,k,:].sum(axis=0),scores[:,k]).statistic for k in range(4)]
        diagnostics.append({'replica':rep,
            'half_kick_by_time':relative_by_time(scores[data['half_kick_sources']],data['half_kick_scores_ps']),
            'half_timestep_by_time':relative_by_time(scores[data['half_timestep_sources']],data['half_timestep_scores_ps']),
            'pair_norm_matrix_relative_asymmetry_by_time':asymmetry,
            'incoming_outgoing_spearman_by_time':incoming_outgoing,
            'no_interactions_max_score_by_time':data['free_scores_ps'].max(axis=0)})
        rank_changes.append({'replica':rep,'force_leader':int(force.argmax()),'force_leader_rank_by_time':ranks,
                             'response_leaders_by_time':scores.argmax(axis=0)})
        label=f'Run {rep+1}'
        ax[0,0].plot(times*1000,correlations,'-o',color=colors[rep],label=label)
        ax[0,1].plot(times*1000,ranks,'-o',color=colors[rep],label=label)
        ax[1,0].scatter(force/force.mean(),scores[:,-1]/scores[:,-1].mean(),s=13,
                        alpha=.5,color=colors[rep],label=label,edgecolor='none')
        ax[1,1].plot(times*1000,scores.std(axis=0)/scores.mean(axis=0)*100,'-o',color=colors[rep],label=label)
    ax[0,0].set(title='A | Predictive value of the original force score',xlabel='Time after impulse (fs)',ylabel='Spearman rank correlation',ylim=(-.3,1.05))
    ax[0,0].axhline(0,color='#777',lw=1,ls=':');ax[0,0].legend(frameon=False)
    ax[0,1].set(title='B | Where the original force leader ranks',xlabel='Time after impulse (fs)',ylabel='Rank among 216 molecules (1 = highest)')
    ax[0,1].invert_yaxis()
    ax[1,0].set(title='C | Initial force score versus response at 200 fs',xlabel='Force score / snapshot mean',ylabel='Finite-time response / snapshot mean')
    ax[1,1].set(title='D | How unequal the finite-time responses remain',xlabel='Time after impulse (fs)',ylabel='Response coefficient of variation (%)')
    for a in ax.flat:
        a.spines[['top','right']].set_visible(False)
    fig.savefig(R/'finite_time_results.png',dpi=180)
    fig.savefig(R/'finite_time_results.svg')
    plt.close(fig)
    result={'times_ps':times,'diagnostics':diagnostics,'rank_changes':rank_changes,
        'scope':'Three conditional snapshots; uncertainty is not estimated by treating the 216 molecules as independent replicates.',
        'verification':'Array dimensions, finite values, nonnegative responses, excluded-source diagonals, and score reconstruction passed.'}
    (D/'analysis.json').write_text(json.dumps(result,indent=2,default=native),encoding='utf-8')
    print(json.dumps(result,indent=2,default=native))


if __name__=='__main__':main()
