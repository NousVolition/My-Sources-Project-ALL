"""Compare midpoint staircase currents to smooth triangular ramps at equal speed."""
import numpy as np
from circuit import to_state
from support import ROOT,P,save,table,interval,evolve

def trajectory(s,leg,bins,low,high,alpha,step=None,max_step=.5):
    t=np.linspace(0,2*leg,2*bins+1);bias=low+(high-low)*(1-abs(t/leg-1))
    if step is None:
        f=lambda u:low+(high-low)*(1-abs(u/leg-1))
        # Explicit split at the slope discontinuity.
        first=evolve(s,f,alpha,leg,t[:bins+1],max_step=max_step)
        second=evolve(first[-1],lambda u:f(u+leg),alpha,leg,t[:bins+1],max_step=max_step)
        return t,bias,np.concatenate([first,second[1:]])
    count=round((high-low)/step);dur=leg/count;states=[s.copy()];current=s.copy()
    for j in range(2*count):
        level=low+step*(j+.5 if j<count else 2*count-j-.5)
        inside=t[(t>j*dur+1e-8)&(t<=(j+1)*dur+1e-8)]-j*dur
        sol=evolve(current,level,alpha,dur,max_step=max_step)
        if len(inside):states.extend(sol.sol(inside).T.reshape((-1,)+s.shape))
        current=sol.y[:,-1].reshape(s.shape)
    return t,bias,np.array(states)

def main():
    c=P['ramps'];rng=np.random.default_rng(c['seed']);phi=rng.uniform(-np.pi,np.pi,(c['starts'],2));s,sign=to_state(phi)
    s=evolve(s,c['low'],c['alpha'],c['burn'],[c['burn']])[-1];rows=[];data={'initial_phase':phi,'burned_state':s};checks=[]
    for leg in c['leg_durations']:
        t,b,ref=trajectory(s,leg,c['bins_per_leg'],c['low'],c['high'],c['alpha'])
        _,_,fine=trajectory(s,leg,c['bins_per_leg'],c['low'],c['high'],c['alpha'],max_step=.25)
        v=2*np.diff(ref[:,:,0],axis=0)/np.diff(t)[:,None];vf=2*np.diff(fine[:,:,0],axis=0)/np.diff(t)[:,None]
        checks.append(float(abs(v-vf).max()));data.update({f'time_{leg}':t,f'bias_{leg}':b,f'smooth_{leg}':ref,f'smooth_fine_{leg}':fine})
        for step in c['current_steps']:
            _,_,q=trajectory(s,leg,c['bins_per_leg'],c['low'],c['high'],c['alpha'],step)
            qv=2*np.diff(q[:,:,0],axis=0)/np.diff(t)[:,None];ci=interval(np.sqrt(np.mean((qv-v)**2,axis=0)))
            rows.append({'leg_duration':leg,'current_increment':step,'mean_run_rms_vs_smooth':ci['mean'],'ci_low':ci['ci95'][0],'ci_high':ci['ci95'][1]});data[f'stair_{leg}_{step}']=q
    np.savez_compressed(ROOT/'data/ramps.npz',**data);table('ramps.csv',rows)
    save('results/ramps.json',{'independent_starts':c['starts'],'regimes':rows,'max_smooth_step_refinement_voltage_error':max(checks),'step_gate_2e_6':max(checks)<2e-6,'design':'Same initial states, current endpoints, elapsed leg duration and voltage bins; only staircase resolution changes.'})
    print('ramps',rows,'refinement',checks,flush=True)
def refine():
    c=P['ramps'];old=np.load(ROOT/'data/ramps.npz');s=old['burned_state'];rows=[];data={}
    for leg in c['leg_durations']:
        t=old[f'time_{leg}'];ref=old[f'smooth_fine_{leg}'];v=2*np.diff(ref[:,:,0],axis=0)/np.diff(t)[:,None]
        for step in (.00625,.003125,.0015625):
            _,_,q=trajectory(s,leg,c['bins_per_leg'],c['low'],c['high'],c['alpha'],step)
            qv=2*np.diff(q[:,:,0],axis=0)/np.diff(t)[:,None];ci=interval(np.sqrt(np.mean((qv-v)**2,axis=0)))
            rows.append({'leg_duration':leg,'current_increment':step,'mean_run_rms_vs_smooth':ci['mean'],'ci_low':ci['ci95'][0],'ci_high':ci['ci95'][1]});data[f'stair_{leg}_{step}']=q
    table('ramps_refined.csv',rows);np.savez_compressed(ROOT/'data/ramps_refined.npz',**data)
    save('results/ramps_refined.json',{'reason':'Post-check numerical follow-up: the original long-leg staircase errors were nonmonotone. Original results are retained. No forecast models were changed.','regimes':rows})
    print('refined ramps',rows,flush=True)
if __name__=='__main__':
    import sys
    refine() if '--refine' in sys.argv else main()
