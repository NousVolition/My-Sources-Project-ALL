"""Audit saved trajectories and endpoints; generate finite-time response figures."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from coupled_model import CoupledFlow
from run_coupled import case_id, write, ROOT, digest

def embedded(h,n=64):
    old=h.shape[-1];m=np.rint(np.fft.fftfreq(old)*old).astype(int)
    out=np.zeros(h.shape[:-2]+(n,n),dtype=complex)
    out[...,m[:,None]%n,m[None,:]%n]=h/old**2
    return out

def norm(h):return np.sqrt(np.sum(abs(h)**2,axis=(-3,-2,-1)))
def delta(h):return h[:,1,:2]-h[:,0,:2]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'output');args=ap.parse_args();out=args.out
    rows=json.loads((out/'results.json').read_text(encoding='utf-8'));p=json.loads((ROOT/'protocol.json').read_text(encoding='utf-8'))
    prov=json.loads((out/'provenance.json').read_text(encoding='utf-8'))
    assert all(digest(ROOT/name)==sha for name,sha in prov['source_sha256'].items())
    hist={};ends=np.load(out/'final-fields.npz');checks=[];scale=rows[0]['initial_impulse_L2']
    max_budget=max_mass=max_div=0.;min_phi=1.;max_phi=0.
    for r in rows:
        key=r['case'];file=out/'trajectories'/(key+'.npz');assert digest(file)==r['trajectory_sha256']
        h=np.load(file)['h'];assert h.shape==(51,2,3,r['n'],r['n']) and np.isfinite(h).all()
        np.testing.assert_array_equal(h[-1],ends[key]);hist[key]=h
        m=CoupledFlow(r['n'],feedback=r['feedback']);e0=[m.diagnostics(h[0,j])['energy'] for j in [0,1]]
        assert abs(r['initial_impulse_L2']/scale-1)<1e-12
        for i,t in enumerate(r['rows']):
            assert abs(t['t']-i*.02)<1e-14
            response=float(norm(embedded(delta(h)[i:i+1]))[0]/scale)
            assert abs(response-t['disturbance_velocity_L2_relative_to_initial_impulse'])<1e-12
            for j,label in enumerate(['baseline','disturbed']):
                d=m.diagnostics(h[i,j]);saved=t[label]
                assert all(abs(d[k]-saved[k])<1e-11 for k in d)
                max_budget=max(max_budget,abs(saved['energy_budget_relative']))
                max_mass=max(max_mass,abs(d['mass']-.02));max_div=max(max_div,d['divergence_max'])
                min_phi=min(min_phi,d['fraction_min']);max_phi=max(max_phi,d['fraction_max'])
        # Independent, lower-order trapezoid budget sanity check from saved times.
        trap=[]
        for j,label in enumerate(['baseline','disturbed']):
            ds=np.array([a[label]['dissipation'] for a in r['rows']]);ts=np.array([a['t'] for a in r['rows']])
            trap.append(float((r['rows'][-1][label]['energy']-e0[j]+np.trapezoid(ds,ts))/e0[j]))
        checks.append({'case':key,'saved_fields':102,'final_field_sha256':hashlib.sha256(ends[key].tobytes()).hexdigest(),'trapezoid_budget_relative':trap})
    def get(n,dt,s,o,f):return hist[case_id(n,dt,s,o,f)]
    def contrast(n,dt,s,f):
        a=get(n,dt,s,'x',f);b=get(n,dt,s,'y',f)
        c=norm(embedded(delta(a))-embedded(delta(b)))/scale
        return {'n':n,'dt':dt,'seed':s,'feedback':f,'rms':float(np.sqrt(np.mean(c*c))),'final':float(c[-1]),'curve':c.tolist()}
    contrasts=[contrast(48,.004,s,f) for s in p['primary_seeds'] for f in [True,False]]
    contrasts += [contrast(n,dt,100,f) for n,dt in p['additional_refinements'] for f in [True,False]]
    refinements=[]
    for (n0,d0),(n1,d1) in [((32,.004),(48,.004)),((48,.004),(64,.004)),((64,.004),(64,.002))]:
        for o in ['x','y']:
            a=get(n0,d0,100,o,True);b=get(n1,d1,100,o,True)
            e=norm(embedded(delta(a))-embedded(delta(b)))/scale
            refinements.append({'from':[n0,d0],'to':[n1,d1],'organization':o,'rms_error':float(np.sqrt(np.mean(e*e))),'max_error':float(e.max()),'initial_error':float(e[0])})
    null=max(c['rms'] for c in contrasts if not c['feedback'])
    uniform=float(norm(embedded(delta(get(48,.004,100,'uniform',True)))-embedded(delta(get(48,.004,100,'uniform',False)))).max()/scale)
    ref=contrast(64,.002,100,True)['rms'];coarse=contrast(48,.004,100,True)['rms']
    err=max(c['rms_error'] for c in refinements if c['from'][0]>=48)
    gates={'all_30_pairs':len(rows)==30,'null':null<1e-11,'uniform':uniform<1e-11,
           'mass':max_mass<1e-11,'divergence':max_div<1e-10,'fraction_bounds':min_phi>=-1e-9 and max_phi<=.03500001,
           'energy_budget':max_budget<1e-5,'initial_fields':max(c['initial_error'] for c in refinements)<1e-12,
           'effect_above_error':ref>10*err,'contrast_refinement':abs(coarse/ref-1)<.01}
    summary={'pairs':len(rows),'trajectories':2*len(rows),'saved_model_time':1.,'gates':gates,'passed':all(gates.values()),
             'contrasts':contrasts,'refinements':refinements,'effect_to_error_ratio':ref/err,
             'relative_contrast_refinement_change':abs(coarse/ref-1),'feedback_off_null':null,'uniform_null':uniform,
             'max_relative_energy_budget':max_budget,'max_mass_error':max_mass,'max_divergence':max_div,
             'fraction_range':[min_phi,max_phi],'endpoint_audits':checks}
    write(out/'analysis.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=np.arange(51)*.02
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for c in contrasts[:8]:
        if c['feedback']:axs[0].plot(t,100*np.array(c['curve']),label=f"Phase seed {c['seed']}")
    axs[0].plot(t,np.zeros(51),'k--',label='Feedback off: zero')
    axs[0].set(xlabel='Model time',ylabel='Organization contrast / initial impulse (%)',title='Disturbance-specific response');axs[0].legend(fontsize=8)
    for r in rows[:4]:axs[1].plot(t,[a['disturbance_velocity_L2_relative_to_initial_impulse'] for a in r['rows']],label=r['organization']+' / '+('feedback on' if r['feedback'] else 'feedback off'))
    axs[1].set(xlabel='Model time',ylabel='Response amplitude / initial impulse',title='Seed 100: decay after the same impulse');axs[1].legend(fontsize=8)
    fig.suptitle('Transported material changes viscous stress under an assumed constitutive law',fontsize=12)
    fig.savefig(out/'response.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for o in ['x','y']:
        c=[r for r in refinements if r['organization']==o]
        axs[0].semilogy(range(3),[max(r['rms_error'],1e-16) for r in c],'o-',label=o)
    axs[0].set_xticks(range(3),['32 to 48','48 to 64','64 half-step'],rotation=20)
    axs[0].axhline(ref/10,color='k',ls='--',label='Effect / 10');axs[0].set(ylabel='RMS field error / initial impulse',title='Seed 100 refinement');axs[0].legend(fontsize=8)
    for ax,o in zip(axs[1:],['x','y']):
        h=get(64,.002,100,o,True);phi=CoupledFlow(64).physical(h[-1,1])[2]
        im=ax.imshow(phi.T,origin='lower',extent=[0,1,0,1],vmin=.005,vmax=.035,cmap='viridis');ax.set(title=f'{o} organization at t=1',xlabel='x',ylabel='y')
    fig.colorbar(im,ax=axs[1:],label='Material volume fraction',shrink=.8)
    fig.savefig(out/'controls.png',dpi=150);plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['contrasts','endpoint_audits']},indent=2))
    if not summary['passed']:raise SystemExit('A protocol gate failed; preserve outputs')

if __name__=='__main__':main()
