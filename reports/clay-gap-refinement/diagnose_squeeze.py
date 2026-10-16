"""Quantify active contact terms in the preserved original peak states."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,numpy as np
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'baseline-fresh'))
from mesh_hug import Mesh,Parameters
from contact_hug import barrier
rows=[]
for n in [144,168]:
 p=R/'baseline-data'/f'shared{n}-one.json';r=json.loads(p.read_text());z=np.load(p.with_suffix('.npz'));m=Mesh(Parameters(**r['parameters']))
 for target in [0,.1,.12,.125,.13,.135,.15,.2]:
  i=int(np.argmin(abs(z['time']-target)));x,rest,l=m.unpack(z['state'][i]);g=m.self_gap(x,False);active=g<m.p.self_width;same=m.si//n==m.sa//n;delta=np.minimum((m.si-m.sa)%n,(m.sa-m.si)%n);e=barrier(g,m.p.self_width,m.p.self_strength)[0]*m.sw
  rows.append(dict(n=n,time=float(z['time'][i]),gap=r['diagnostics'][i]['gap_arclength_rms'],original_gap=r['diagnostics'][i]['surface_gap_rms'],self_energy=float(sum(e)),self_active=int(sum(active)),same_surface=int(sum(active&same)),within_three_reference_angle_steps=int(sum(active&same&(delta<4))),minimum_self=float(min(g))))
(R/'data/squeeze-diagnostic.json').write_text(json.dumps(dict(rows=rows,interpretation='Localization/contact diagnostic only; not an isolated causal attribution'),indent=2)+'\n');print(json.dumps(dict(states=len(rows))))
