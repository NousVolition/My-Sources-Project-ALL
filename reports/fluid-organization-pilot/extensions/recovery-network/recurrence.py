"""Descriptive replay: returns near the complete pre-shutoff velocity field.

Added after initial network analysis, before fluid prediction analysis. No model
selection uses these outcomes. Replay energies must match original recordings.
"""
import numpy as np
from core import PilotFlow,initial,energy,L
from recovery import rk4
from common import ROOT,P,save

def return_events(t,distance,exit_radius=.2,return_radius=.1,residence=.2):
    """Exit/return hysteresis; return must persist through a full sampled window."""
    events=[];outside=False
    for i in range(len(t)):
        if not outside and distance[i]>=exit_radius:outside=True
        if outside and distance[i]<=return_radius:
            j=np.searchsorted(t,t[i]+residence-1e-9)
            if j<len(t) and np.max(distance[i:j+1])<=return_radius:
                events.append(float(t[i]));outside=False
    return events

def main():
    c=P['fluid'];rows=[];dest=ROOT/'data'/'recurrence';dest.mkdir(parents=True,exist_ok=True)
    for seed in c['test_seeds']:
        recorded=np.load(ROOT/'data'/'fluid'/f'base_{seed}.npz');h0=recorded['initial_h'];f=PilotFlow(c['n'],c['nu'])
        force=initial(f,7788,kind='beltrami')*c['forcing_rms'];hs=[h0.copy() for _ in c['arms']];dt=c['dt'];stride=round(c['save_dt']/dt)
        distances=[];es=[];e0=energy(f,h0)
        for j in range(round(c['end']/dt)+1):
            if j%stride==0:
                distances.append([np.sqrt(energy(f,q-h0)/e0) for q in hs]);es.append([energy(f,q) for q in hs])
            if j<round(c['end']/dt):
                for a,arm in enumerate(c['arms']):hs[a],_,_=rk4(f,hs[a],j*dt,dt,force,arm,c['ramp_time'])
        error=float(np.max(abs(np.array(es)-recorded['energy']))/e0)
        if error>1e-12:raise AssertionError('Replay differs from primary trajectory')
        d=np.array(distances);t=recorded['time'];np.savez_compressed(dest/f'seed{seed}.npz',time=t,distance_to_initial=d)
        for a,arm in enumerate(c['arms']):
            row={'seed':seed,'arm':arm,'relative_l2_final':float(d[-1,a]),'replay_energy_error':error,
                 'exited_radius_0_2':bool(np.max(d[:,a])>=.2),'return_times':return_events(t,d[:,a])}
            row['radius_sensitivity']={str(r):return_events(t,d[:,a],2*r,r) for r in (.05,.1,.2)};rows.append(row)
        print(f'recurrence {seed} verified',flush=True)
    save(ROOT/'results'/'recurrence.json',{'definition':'Relative spatial L2 distance from velocity field at shutoff; exit 0.2, return 0.1, residence 0.2. Threshold sensitivity 0.05/0.1/0.2; not a search over all historical states.',
            'records':rows,'total_returns':sum(len(r['return_times']) for r in rows),'n_independent_test_runs':len(c['test_seeds'])})

if __name__=='__main__':main()
