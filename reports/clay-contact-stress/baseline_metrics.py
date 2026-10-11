"""Save actual molding, reciprocal-coupling and numerical stress results."""
from pathlib import Path
from dataclasses import asdict,replace
import json,hashlib,time
import numpy as np
from scipy.integrate import trapezoid
from clay_hug import *
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def save(n,d):(DATA/n).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def summarize(mesh,t,y,burn=4):
    diagnostics=[mesh.diagnostics(s) for s in y]
    e=np.array([d['material_energy'] for d in diagnostics]);l=np.array([d['lorenz_energy'] for d in diagnostics]);native,visc,plastic,port=y[:,-4:].T
    material_res=e-e[0]+visc+plastic-port
    lorenz_res=l-l[0]-native+(port if mesh.p.feedback else 0)
    total_res=e+l-e[0]-l[0]-native+visc+plastic-(0 if mesh.p.feedback else port)
    def scaled(res,energy,work):return float(max(abs(res))/max(1,max(abs(energy)),max(abs(work))))
    mask=t>=min(burn,t[-1]/2);tt=t[mask];ly=y[mask,-7:-4]
    mean=lambda v:float(trapezoid(v,tt)/(tt[-1]-tt[0])) if len(tt)>1 else float(v[0])
    zmean=mean(ly[:,2]);imprint=np.array([d['imprint'] for d in diagnostics])
    velocity=np.array([mesh.rhs(a,b)[:mesh.position_size] for a,b in zip(t,y)])
    rate=np.array([mesh.rhs(a,b)[-1] for a,b in zip(t,y)])
    return dict(final=diagnostics[-1],initial=diagnostics[0],
        lorenz_rms=float(np.sqrt(mean(np.sum(ly[:,:2]**2,axis=1)))),z_std=float(np.sqrt(max(0,mean((ly[:,2]-zmean)**2)))),
        mean_imprint=mean(imprint[mask]),mean_contact_fraction=mean(np.array([d['contact_fraction'] for d in diagnostics])[mask]),
        mesh_motion_rms=float(np.sqrt(mean(np.sum(velocity[mask]**2,axis=1)/mesh.n))),
        maximum_penetration=max(d['penetration'] for d in diagnostics),maximum_wall_overlap=max(d['wall_overlap'] for d in diagnostics),
        minimum_cell_area_ratio=min(d['min_area_ratio'] for d in diagnostics),
        maximum_area_change=max(abs(d['area_ratio']-1) for d in diagnostics),
        material_budget_error=scaled(material_res,e,port),lorenz_budget_error=scaled(lorenz_res,l,native),
        total_budget_error=scaled(total_res,e+l,native),plastic_dissipation=float(plastic[-1]),viscous_dissipation=float(visc[-1]),
        net_port_work=float(port[-1]),sampled_return_work=float(trapezoid(np.maximum(-rate,0),t)),
        sampled_forward_work=float(trapezoid(np.maximum(rate,0),t)))

def main():
    DATA.mkdir(exist_ok=True)
    assert json.loads((DATA/'validation.json').read_text())['all_passed']
    rows=[]
    def run(name,p,end=16,dt=.001,initial=(1.,1.,1.),group='paired',state=None):
        start=time.perf_counter();m,t,y,f=integrate(p,end=end,dt=dt,initial=initial,state=state)
        row=dict(name=name,group=group,parameters=asdict(p),initial=list(initial),end=end,dt=dt,failure=f,seconds=time.perf_counter()-start,
                 metrics=summarize(m,t,y) if len(t)>1 else None)
        np.savez_compressed(DATA/(name+'.npz'),time=t,state=y,triangles=m.tri,initial_positions=m.initial,outline=object_outline(p))
        rows.append(row);save('runs.json',rows)
        print(name,'completed',t[-1],'failure',f,'contact',row['metrics']['final']['contact_fraction'] if row['metrics'] else None,flush=True)
        return m,t,y,row
    base=Parameters(preload=.35)
    configs=[]
    for shape in ['circle','ellipse','lobed']:
        for seed,initial in enumerate([(1.,1.,1.),(-8.,8.,27.)]):
            for feedback in [False,True]:
                p=replace(base,object_shape=shape,feedback=feedback)
                configs.append((f'{shape}-s{seed}-'+('two' if feedback else 'one'),p,initial))
    for alpha in [.15,1.,2.]:
        for feedback in [False,True]:configs.append((f'alpha{alpha:g}-'+('two' if feedback else 'one'),replace(base,alpha=alpha,feedback=feedback),(1.,1.,1.)))
    save('protocol.json',dict(parameters=asdict(base),paired_cases=len(configs),seeds=[[1,1,1],[-8,8,27]],
        record_step=.04,paired_end=16,statistics_burn=4,dt=.001,
        sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['clay_hug.py','run_stress.py','validate_model.py']},
        scope='Dimensionless 2-D cohesive viscoplastic mesh and prescribed rigid objects. Not calibrated real clay, 3-D material, water or Navier-Stokes flow.',
        port='Force = -(alpha/10)*X*grad(A); reciprocal Xdot += (alpha/10)*Adot / lorenz_capacity; A=sum(reference nodal area*radius)',
        targets=dict(short_state_scaled_error=1e-4,energy_budget=1e-3,radial_overlap=.02,minimum_area_ratio=.1)))
    for name,p,initial in configs:run(name,p,initial=initial)
    # Long finite windows and a changed reservoir capacity test reciprocal influence.
    for feedback in [False,True]:run('long-'+('two' if feedback else 'one'),base if feedback else replace(base,feedback=False),end=40,group='long')
    for capacity in [.0025,.04]:run(f'capacity-{capacity:g}',replace(base,lorenz_capacity=capacity),group='capacity')
    # Check time resolution in the same configuration, on a short useful horizon.
    checks=[]
    for name,p in [('ellipse',base),('lobed',replace(base,object_shape='lobed')),('strong',replace(base,alpha=2))]:
        a=run('step-'+name+'-base',p,end=4,group='convergence')
        b=run('step-'+name+'-fine',p,end=4,dt=.0005,group='convergence')
        scale=np.maximum(1,np.max(abs(b[2][:,:-4]),axis=0));error=float(np.max(abs(a[2][:,:-4]-b[2][:,:-4])/scale))
        checks.append(dict(name=name,scaled_state_error=error,passed=bool(error<1e-4 and a[3]['failure'] is None and b[3]['failure'] is None)))
    # Different adaptive integrator uses the identical continuous equations.
    m,t,y,f=integrate(base,end=4,dt=.003,method='DOP853')
    np.savez_compressed(DATA/'independent.npz',time=t,state=y)
    stored=np.load(DATA/'step-ellipse-fine.npz');scale=np.maximum(1,np.max(abs(y[:,:-4]),axis=0))
    error=float(np.max(abs(y[:,:-4]-stored['state'][:,:-4])/scale))
    checks.append(dict(name='independent-DOP853',scaled_state_error=error,passed=bool(error<1e-4)))
    save('solver-checks.json',checks)
    # Spatial resolution under steady forming avoids chaotic phase drift masking mesh error.
    for angles,layers,dt in [(24,3,.001),(48,5,.0005),(72,7,.00025)]:
        run(f'mesh-{angles}',replace(base,alpha=0,angles=angles,layers=layers),end=12,dt=dt,group='mesh')
    for contact,dt in [(25,.001),(100,.001),(400,.0005)]:
        run(f'contact-{contact}',replace(base,alpha=0,contact=contact),end=12,dt=dt,group='contact')
    # Plastic vs elastic control, followed by release with the object still present.
    for plastic in [False,True]:
        p=replace(base,alpha=0,plastic=plastic)
        a=run('forming-'+('clay' if plastic else 'elastic'),p,end=12,group='memory')
        run('release-'+('clay' if plastic else 'elastic'),replace(p,preload=0),end=30,state=a[2][-1],group='memory')
    # Adversarial time steps and physical parameters; preserve failed attempts.
    for key,p in [('strong',replace(base,alpha=4)),('fast',replace(base,viscosity=.15)),
        ('stiff',replace(base,bulk=50,contact=1000)),('high-rho',replace(base,rho=100,alpha=1))]:
        run('stress-'+key+'-coarse',p,end=8,dt=.008,group='stress')
        run('stress-'+key+'-resolved',p,end=8,dt=.00025,group='stress')
    summary=dict(configurations=len(rows),all_short_solver_checks_passed=all(c['passed'] for c in checks),
        failures=[dict(name=r['name'],failure=r['failure']) for r in rows if r['failure']],
        largest_completed_material_budget_error=max(r['metrics']['material_budget_error'] for r in rows if r['failure'] is None),
        largest_completed_total_budget_error=max(r['metrics']['total_budget_error'] for r in rows if r['failure'] is None),
        maximum_overlap=max(r['metrics']['maximum_penetration'] for r in rows if r['failure'] is None),
        note='All configurations are recorded, including failed coarse steps and parameter cases that miss geometric targets. Finite windows and two starting states do not establish global stability or statistical confidence intervals.')
    save('summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    if not summary['all_short_solver_checks_passed']:raise SystemExit(1)
if __name__=='__main__':main()
