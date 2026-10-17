"""Preserve the executed solver; create a packed streaming wrapper for 128^3.

Only starting-field storage and branch scheduling change. Retained coefficients
remain complex128; n=128 normalization is exact binary scaling. No equations,
preparation, observation formula, precision or physical parameters change.
"""
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    destination=ROOT/'bounded_final_experiment.py'
    if destination.exists():raise RuntimeError('Bounded source already exists; do not overwrite an executed version')
    s=(ROOT/'final_experiment.py').read_text(encoding='utf-8')
    s=s.replace('from corrected_solver import random_field, interventions','from corrected_solver import random_field, interventions\nfrom final_experiment import spinup as original_spinup')
    s=s.replace("'corrected_observation.py','final_experiment.py']","'corrected_observation.py','final_experiment.py','bounded_final_experiment.py']")
    s=s.replace('sr,spin_meta=spinup(seed,regime,backend)','sr,spin_meta=original_spinup(seed,regime,backend)')
    start=s.index('    states=[];labels=[]',s.index('def run_one'))
    stop=s.index('    times=np.arange',start)
    s=s[:start]+'''    states=[];labels=[]
    for i,(arm,hist) in enumerate((a,j) for a in P['arms'] for j in range(2)):
        states.append(s.pack(bases[i:i+1])[0]);labels.append(dict(arm=arm,history=hist,method='unprobed',base=i))
        for j,method in enumerate(P['probe_methods']):
            pulse=xp.einsum('i,icxyz->cxyz',xp.asarray(coeff[i,j]),q)
            states.append(s.pack((bases[i]+pulse)[None])[0]);labels.append(dict(arm=arm,history=hist,method=method,base=i))
    start_fields=np.stack(states);del states,pulse
    batch=6 if n>=128 else (16 if n>=64 else len(labels));rows=[];ends=[]
    if backend=='gpu':xp.get_default_memory_pool().free_all_blocks()
'''+s[stop:]
    s=s.replace('h=initial_branches[offset:offset+batch].copy();b=len(h);','h=s.unpack(start_fields[offset:offset+batch]);b=len(h);')
    s=s.replace('end_fields=np.concatenate(ends);start_fields=s.pack(initial_branches)','end_fields=np.concatenate(ends)')
    assert 'initial_branches' not in s
    destination.write_text(s,encoding='utf-8')
    print('Created separate bounded source; original executed source unchanged.')

if __name__=='__main__':main()
