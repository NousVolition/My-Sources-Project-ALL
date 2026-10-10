"""Independent short-horizon checks at every Lorenz rho in the new sweep."""
from pathlib import Path
import json
import numpy as np
from lorenz_hug import Drive,model,integrate
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'lorenz-data'
def main():
    rows=json.loads((DATA/'screen.json').read_text());checks=[]
    for row in rows:
        if row['hug']['damping']!=.05 or row['hug']['memory_coupling']!=1.5 or row['drive']['amplitude']!=.65:continue
        d=Drive(**row['drive']);p=model.Parameters(**row['hug']);a=np.load(DATA/(row['name']+'.npz'))
        t=a['time'][:501];base=a['state'][:501];start=base[0,:3]
        sol,opened,_=integrate(d,p,end=5,initial=start,method='Radau',rtol=1e-11,max_step=.005)
        ref=sol.sol(t).T;scale=np.maximum(1,np.max(abs(ref[:,:6]),axis=0));error=float(np.max(abs(base[:,:6]-ref[:,:6])/scale))
        checks.append(dict(name=row['name'],rho=d.rho,scaled_error=error,target=1e-5,passed=error<1e-5))
        np.savez_compressed(DATA/(row['name']+'-independent.npz'),time=t,state=ref)
    (DATA/'sweep-independent.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
    # Exercise an actual upward opening crossing after startup, not only an initially open state.
    d=Drive(mean=.75,amplitude=.65);p=model.Parameters(r=-.25,memory_coupling=1.5);events=[];arrays={}
    for method,step in [('DOP853',.02),('Radau',.02),('DOP853',.01)]:
        sol,opened,_=integrate(d,p,end=2,initial=[1,1,1],method=method,max_step=step)
        events.append(opened);arrays[method+str(step)]=sol.sol(np.linspace(0,2,401)).T
    assert all(t is not None and 0<t<2 for t in events)
    event=dict(opening_times=events,maximum_difference=max(events)-min(events),target=1e-7,passed=max(events)-min(events)<1e-7)
    (DATA/'event-check.json').write_text(json.dumps(event,indent=2)+'\n');np.savez_compressed(DATA/'event-check.npz',time=np.linspace(0,2,401),**arrays)
    print(json.dumps(event))
    if not all(x['passed'] for x in checks) or not event['passed']:raise SystemExit(1)
if __name__=='__main__':main()
