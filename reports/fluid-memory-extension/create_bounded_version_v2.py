"""Build the final memory-bounded 128^3 wrapper from the unchanged solver.

Release unneeded preparation arrays once all exact starting coefficients and
the probe basis are on the CPU. Eight branches then fit without retaining the
entire 56-branch Fourier grid on the GPU. No mathematical formula changes.
"""
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    destination=ROOT/'bounded_final_experiment_v2.py'
    if destination.exists():raise RuntimeError('Bounded v2 source exists; do not overwrite an executed version')
    s=(ROOT/'final_experiment.py').read_text(encoding='utf-8')
    s=s.replace('from corrected_solver import random_field, interventions','from corrected_solver import random_field, interventions\nfrom final_experiment import spinup as original_spinup')
    s=s.replace("'corrected_observation.py','final_experiment.py']","'corrected_observation.py','final_experiment.py','bounded_final_experiment_v2.py']")
    s=s.replace('sr,spin_meta=spinup(seed,regime,backend)','sr,spin_meta=original_spinup(seed,regime,backend)')
    start=s.index('    states=[];labels=[]',s.index('def run_one'))
    stop=s.index('    times=np.arange',start)
    s=s[:start]+'''    states=[];labels=[]
    for i,(arm,hist) in enumerate((a,j) for a in P['arms'] for j in range(2)):
        states.append(s.pack(bases[i:i+1])[0]);labels.append(dict(arm=arm,history=hist,method='unprobed',base=i))
        for j,method in enumerate(P['probe_methods']):
            pulse=xp.einsum('i,icxyz->cxyz',xp.asarray(coeff[i,j]),q)
            states.append(s.pack((bases[i]+pulse)[None])[0]);labels.append(dict(arm=arm,history=hist,method=method,base=i))
    start_fields=np.stack(states);probe_basis=s.pack(q)
    del states,pulse,variants,bases,scr,shifted,initial,writer,q,force,ff,h
    batch=8 if n>=128 else (16 if n>=64 else len(labels));rows=[];ends=[]
    if backend=='gpu':xp.get_default_memory_pool().free_all_blocks()
'''+s[stop:]
    s=s.replace('h=initial_branches[offset:offset+batch].copy();b=len(h);','h=s.unpack(start_fields[offset:offset+batch]);b=len(h);')
    s=s.replace('end_fields=np.concatenate(ends);start_fields=s.pack(initial_branches)','end_fields=np.concatenate(ends)')
    s=s.replace('probe_basis=s.pack(q),pulse_coefficients=coeff','probe_basis=probe_basis,pulse_coefficients=coeff')
    assert 'initial_branches' not in s
    destination.write_text(s,encoding='utf-8')
    print('Created v2 packed wrapper with preparation arrays released before readout.')

if __name__=='__main__':main()
