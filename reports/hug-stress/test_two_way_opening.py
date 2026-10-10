"""Threshold and opening-release checks for the proposed two-way closure."""
from pathlib import Path
from dataclasses import asdict
import hashlib,json
import numpy as np
import two_way_hug as tw
from lorenz_hug import Drive,pressure,model

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'two-way-data'
def main():
    rows=[];p=model.Parameters(r=-.25,memory_coupling=1.5)
    starts=[np.load(ROOT/'lorenz-data/chaotic.npz')['state'][0,:3],np.array([1.,1.,27.])]
    for seed,start in enumerate(starts):
        for initial_pressure in [.99,1.,1.01]:
            mean=initial_pressure-.6*np.tanh((start[2]-27)/10)
            d=Drive(mean=mean,amplitude=.6)
            for g in [0.,8.]:
                a=tw.integrate(d,p,g,end=8,initial=start)
                b=tw.integrate(d,p,g,end=8,initial=start,rtol=1e-11,max_step=.005)
                t=np.linspace(0,8,1601);y=a.sol(t).T;z=b.sol(t).T
                scale=np.maximum(1,np.max(abs(z[:,:6]),axis=0))
                error=float(np.max(abs(y[:,:6]-z[:,:6])/scale))
                same=(a.opening is None)==(b.opening is None)
                event_error=abs(a.opening-b.opening) if a.opening is not None and b.opening is not None else 0.
                eta=np.array([tw.resistance(tt,mm,g,a.opening) for tt,mm in zip(t,y[:,5])])
                after=t>=a.opening+.8 if a.opening is not None else np.zeros(len(t),bool)
                release=float(max(eta[after])) if np.any(after) else None
                row=dict(name=f'opening-{len(rows):02}',seed=seed,initial=start.tolist(),drive=asdict(d),strength=g,
                    requested_initial_pressure=initial_pressure,actual_initial_pressure=float(pressure(start[2],d)),
                    opening=a.opening,refined_opening=b.opening,scaled_error=error,event_error=event_error,
                    resistance_after_release=release,peak_imprint=float(max(y[:,5])),final_imprint=float(y[-1,5]),
                    feedback_loss=float(y[-1,8]),passed=bool(same and error<1e-5 and event_error<1e-6 and release==0.))
                rows.append(row)
                np.savez_compressed(DATA/(row['name']+'.npz'),time=t,state=y,refined=z,pressure=pressure(y[:,2],d),resistance=eta)
    result=dict(cases=rows,all_passed=all(r['passed'] for r in rows),
        maximum_error=max(r['scaled_error'] for r in rows),max_event_error=max(r['event_error'] for r in rows),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Opening is latched: resistance vanishes after 0.8 model time and never returns. Imprint persists and pressure continues. This tests the proposed rule, not physical venting.')
    (DATA/'opening-tests.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
    if not result['all_passed']:raise SystemExit(1)

if __name__=='__main__':main()
