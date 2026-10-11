"""Prespecified secondary full-field response projection on the writer template."""
import numpy as np
from run import ROOT,P,dump
from analyze import load_run,branch_index,interval


def main():
    rows=[]
    for seed in P['test_seeds']:
        a,m=load_run(ROOT/'data',seed)
        weight=np.where(a['modes'][:,2]==0,1.,2.);writer=a['writer']
        # Writer stores one vector-field template, normalized Fourier coefficients.
        writer=writer.reshape(3,-1);wn=np.sqrt(np.sum(abs(writer)**2*weight))
        for t in P['horizons']:
            ri=np.flatnonzero(np.isclose(a['response_time'],t))[0]
            for arm in ['retained']+['scramble'+str(i) for i in range(P['phase_scramble_replicates'])]+['sham','reset']:
                response=np.stack([a['response_fields'][ri,branch_index(m,arm,h,0)//2] for h in range(2)])
                difference=response[0]-response[1]
                projection=float(np.sum((difference*np.conj(writer)).real*weight)/wn/P['probe_velocity_rms'])
                norm=float(np.sqrt(np.sum(abs(difference)**2*weight)))
                cosine=float(np.sum((difference*np.conj(writer)).real*weight)/max(wn*norm,1e-30))
                rows.append(dict(seed=seed,horizon=t,arm=arm,writer_projection_per_probe_rms=projection,writer_cosine=cosine))
    summary={}
    for t in P['horizons']:
        summary[str(t)]={}
        for arm in ['retained','sham','reset']+['scramble'+str(i) for i in range(4)]:
            v=[r['writer_projection_per_probe_rms'] for r in rows if r['horizon']==t and r['arm']==arm]
            summary[str(t)][arm]=interval(v)
    dump(ROOT/'results/template-projection-results.json',dict(horizons=summary,per_flow=rows,interpretation='Exploratory Eulerian alignment with the fixed writing-force template. A signed projection is not a causal leader or a unique state-memory carrier.'))
    import json
    causal=json.loads((ROOT/'results/causal-results.json').read_text());field={}
    for t in P['horizons']:
        r=[v for v in causal['per_flow'] if v['horizon']==t]
        field[str(t)]={}
        for arm in ['retained','sham','reset']:
            field[str(t)][arm]=interval([v['field_response'][arm]['history_difference'] for v in r])
        scramble=[np.mean([v['field_response']['scramble'+str(i)]['history_difference'] for i in range(4)]) for v in r]
        field[str(t)]['scrambled']=interval(scramble)
        field[str(t)]['retained_minus_scrambled']=interval(np.array(field[str(t)]['retained']['per_flow'])-scramble)
    dump(ROOT/'results/response-field-summary.json',dict(horizons=field,definition='Volume-RMS norm of the difference between the two histories\' probe-minus-unprobed velocity fields, divided by probe RMS amplitude. Includes all retained Fourier modes, not only observer points. Exploratory secondary outcome.'))


if __name__=='__main__':main()
