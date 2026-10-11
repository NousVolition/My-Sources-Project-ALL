"""Audit work done by the identical velocity impulse in each present flow."""
import numpy as np
from run import ROOT,P,dump
from analyze import load_run,branch_index,interval

if __name__=='__main__':
    rows=[]
    for seed in P['test_seeds']:
        a,m=load_run(ROOT/'data',seed,keys=['start_fields','modes','probe'])
        w=np.where(a['modes'][:,2]==0,1.,2.)
        for arm in dict.fromkeys(b['arm'] for b in m['branches']):
            energy=[]
            for history in [0,1]:
                h=a['start_fields'][branch_index(m,arm,history,0)];p=a['probe']
                # Stable exact kinetic-energy increment <u,p> + ||p||^2 / 2.
                energy.append(float(np.sum((h*np.conj(p)).real*w)+.5*np.sum(abs(p)**2*w)))
            rows.append(dict(seed=seed,arm=arm,energy_injection=energy,history_energy_injection_difference=energy[0]-energy[1]))
    summary={arm:interval([r['history_energy_injection_difference'] for r in rows if r['arm']==arm]) for arm in dict.fromkeys(r['arm'] for r in rows)}
    dump(ROOT/'results/impulse-energy-audit.json',dict(per_flow=rows,summary=summary,
        interpretation='The fixed impulse has identical velocity shape/amplitude, not identical work. Work depends on present velocity alignment, which is one possible phase-sensitive response pathway. This audit does not isolate that pathway or establish mediation; equal-work and alternate-location probes remain unrun.'))
