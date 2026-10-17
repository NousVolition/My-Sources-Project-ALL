"""Descriptive endpoint velocity-response footprint; not a collective-memory test.

Evaluate saved Fourier polynomials on the same 128^3 quadrature grid. This is
postprocessing only, with no new physics, trajectory integration or fitted cutoff.
"""
import numpy as np
from scipy import fft
from final_experiment import ROOT,P,dump,digest
from analysis import read,effective_regimes,interval
from followup_analysis import path,common_fields
GRID=128
RADII=np.array([.5,1.,1.5,2.,2.5,3.])
def footprint(coeff,modes,center):
    axes=np.arange(GRID)*(2*np.pi/GRID)
    coordinates=np.meshgrid(axes,axes,axes,indexing='ij',sparse=True)
    radius2=sum(((coordinates[j]-center[j]+np.pi)%(2*np.pi)-np.pi)**2 for j in range(3))
    index=tuple(modes[:,j]%GRID for j in range(3))
    rows=[]
    for c in coeff:
        h=np.zeros((3,GRID,GRID,GRID//2+1),complex)
        h[:,index[0],index[1],index[2]]=c*GRID**3
        u=fft.irfftn(h,s=(GRID,)*3,axes=(-3,-2,-1),workers=2)
        energy=np.sum(u*u,axis=0);total=float(energy.sum())
        weights=np.where(modes[:,2]==0,1.,2.);spectral=float(np.sum(abs(c)**2*weights))
        if abs(total/GRID**3-spectral)>1e-12*max(spectral,1e-20):raise RuntimeError('Spatial Parseval consistency failed')
        if spectral<1e-28:
            rows.append(dict(rms=np.sqrt(spectral),radial_fraction=None,outside_radius2=None,rms_periodic_radius=None));continue
        fractions=[float(energy[radius2<r*r].sum()/total) for r in RADII]
        rows.append(dict(rms=np.sqrt(spectral),radial_fraction=fractions,outside_radius2=1-fractions[3],rms_periodic_radius=float(np.sqrt(np.sum(energy*radius2)/total))))
    return rows
def assess():
    rows=[]
    for regime,cfg in effective_regimes().items():
        seeds=P['smooth_seeds'] if regime=='smooth' else P['forced_seeds']
        for layout in P['layouts'] if regime=='smooth' else ['coincident']:
            for seed in seeds:
                file=path('common-probe-data',seed,regime,cfg['grids'][-1],cfg['dt'],layout)
                a,m=read(file,['modes','pulse','end_fields','base_end_fields'])
                response=common_fields(a).reshape(7,2,3,-1);contrast=response[:,0]-response[:,1]
                center=np.array([np.pi]*3)+P['layouts'][layout]['probe_shift']
                values=footprint(np.concatenate([a['pulse'][None],contrast]),a['modes'],center)
                rows.append(dict(seed=seed,regime=regime,layout=layout,forcing='on' if cfg['injection'] else 'none',input=values[0],arms=dict(zip(P['arms'],values[1:])),source_raw_sha256=m['raw_sha256']))
        if regime!='smooth':
            for seed in seeds[:2]:
                file=path('common-unforced-data',seed,regime,cfg['grids'][-1],cfg['dt'])
                a,m=read(file,['modes','pulse','end_fields']);response=common_fields(a).reshape(7,2,3,-1)
                values=footprint(np.concatenate([a['pulse'][None],response[:,0]-response[:,1]]),a['modes'],[np.pi]*3)
                rows.append(dict(seed=seed,regime=regime,layout='coincident',forcing='off',input=values[0],arms=dict(zip(P['arms'],values[1:])),source_raw_sha256=m['raw_sha256']))
    groups={}
    for key in sorted({(r['regime'],r['layout'],r['forcing']) for r in rows}):
        rr=[r for r in rows if (r['regime'],r['layout'],r['forcing'])==key]
        group={}
        for arm in ['input','retained','scrambled','matched_translation']:
            values=[]
            for r in rr:
                if arm=='input':v=r['input']
                elif arm=='scrambled':
                    scr=[r['arms']['scramble'+str(j)] for j in range(4)]
                    v={metric:np.mean([s[metric] for s in scr],axis=0) for metric in ['rms','outside_radius2','rms_periodic_radius','radial_fraction']}
                else:v=r['arms'][arm]
                values.append(v)
            group[arm]={metric:interval([v[metric] for v in values]) for metric in ['rms','outside_radius2','rms_periodic_radius']}
            group[arm]['radial_mean']=np.mean([v['radial_fraction'] for v in values],axis=0).tolist()
        groups['/'.join(key)]=group
    result=dict(description='Eulerian footprint of the difference between histories of probed-minus-unprobed endpoint velocity, relative to initial probe position',evaluation_grid=GRID,radii=RADII.tolist(),time=P['horizon'],groups=groups,per_flow=rows,
        limitations='Endpoint only. Input tails, advection, pressure projection and sustaining-forcing feedback can broaden the footprint. Nonlocal velocity response is not proof of persistent collective memory or a leader. Fractions use a fixed quadrature grid, not a new NS integration.',
        sources={f:digest(ROOT/f) for f in ['spatial_assessment.py','spatial-analysis-plan.json']})
    dump(ROOT/'results'/'spatial-analysis.json',result);print('Descriptive spatial assessment: '+str(len(rows))+' flow/layout conditions',flush=True);return result
if __name__=='__main__':assess()
