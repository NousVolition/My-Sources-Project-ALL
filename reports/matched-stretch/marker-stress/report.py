"""Rebuild the marker stress report and figures from the saved tracking records."""
from pathlib import Path
import os,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
L=6.
def wrap(x):return (x+3)%6-3
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    data=json.loads((ROOT/'results.json').read_text());runs=data['runs'];raw={n:np.load(ROOT/(n+'.npz'),allow_pickle=False) for n in runs};threshold=.05*L/64
    for n,z in raw.items():
        assert sha(ROOT/(n+'.npz'))==runs[n]['trajectory_sha256']
        assert np.isfinite(z['positions']).all() and np.isfinite(z['velocities']).all()
        assert np.array_equal(z['positions'][:,500:],z['positions'][:,:5]) and z['time'][-1]==.3
        assert np.allclose(z['positions'][0,:250],data['initial_positions'],rtol=0,atol=0)
        assert np.allclose(np.linalg.norm(wrap(z['positions'][0,250:500]-z['positions'][0,:250]),axis=1),.0009375,rtol=1e-12)
    comparisons={}
    specs=[('64 tracking step: coarse to medium','tracking_coarse','tracking_medium'),('64 tracking step: medium to fine','tracking_medium','reference'),('128 tracking step: medium to fine','grid128_medium','grid128'),('128 tracking step: fine to refined','grid128','grid128_refined'),('Fluid timestep halved','fluid_half','reference'),('Velocity snapshots coarsened','saved_time_coarse','reference'),('64 to 128 grid','reference','grid128_refined')]
    for label,a,b in specs:
        za,zb=raw[a],raw[b];ix=np.array([np.argmin(abs(zb['time']-t)) for t in za['time']]);assert np.max(abs(zb['time'][ix]-za['time']))<1e-12
        delta=wrap(za['positions'][:,:250]-zb['positions'][ix,:250]);d=np.linalg.norm(delta,axis=2);early=za['time']<=.1
        drift=za['time'][:,None,None]*(np.array(runs[a]['mean_fluid_velocity'])-runs[b]['mean_fluid_velocity']);adjusted=np.linalg.norm(wrap(za['positions'][:,:250]-zb['positions'][ix,:250]-drift),axis=2)
        comparisons[label]={'a':a,'b':b,'times':za['time'].tolist(),'maximum_by_time':d.max(axis=1).tolist(),'rms_by_time':np.sqrt(np.mean(d*d,axis=1)).tolist(),'maximum_through_0_1':float(d[early].max()),'maximum_through_0_3':float(d.max()),'final_rms':float(np.sqrt(np.mean(d[-1]**2))),'screen_through_0_1':bool(d[early].max()<=threshold),'screen_through_0_3':bool(d.max()<=threshold),'mean_drift_adjusted_maximum':float(adjusted.max()),'groups':{g:{'final_rms':float(np.sqrt(np.mean(d[-1,j*125:(j+1)*125]**2))),'maximum':float(d[:,j*125:(j+1)*125].max())} for j,g in enumerate(data['groups'])}}
    summary={'scope':'505 passive tracks in each of 8 reconstructions through 0.30; 250 main labels,250 slightly shifted companions,5 identical clones; no new fluid evolution','path_screen_threshold':threshold,'comparisons':comparisons,'endpoints':{n:{g:runs[n]['metrics'][-1][j] for j,g in enumerate(data['groups'])} for n in ['reference','grid128_refined']},'initials':{n:{g:runs[n]['metrics'][0][j] for j,g in enumerate(data['groups'])} for n in ['reference','grid128_refined']}}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    colors=['#246bce','#da7130'];titles=['Central neighborhood','Periodic-join neighborhood']
    fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for g in (0,1):
        for n,style,color,label in [('reference','-',colors[g],'64 grid'),('grid128_refined','--','#252525','128 grid'),('fluid_half',':','#2d9960','64 fluid half step')]:
            r=runs[n];t=r['times'];m=r['metrics']
            vals=[[v[g]['neighbor_distance_quantiles'][2] for v in m],[v[g]['velocity_deviation_rms'] for v in m],[v[g]['initial_neighbor_retention']*100 for v in m]]
            for ax,y in zip(axs[g],vals):ax.plot(t,y,style,color=color,label=label)
        for ax,title,ylabel in zip(axs[g],['Neighbor spacing','Velocity spread','Neighbors retained'],['Median distance','RMS velocity about cloud mean','Percent']):ax.set(title=titles[g]+'\n'+title,xlabel='Model time',ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Individual fluid markers under stress | original start | grid-sensitive paths');fig.savefig(ROOT/'individual-motion.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for label,c in comparisons.items():
        ax=axs[0] if 'tracking' in label else axs[1];ax.semilogy(c['times'][1:],np.maximum(c['maximum_by_time'][1:],1e-12),label=label)
    for ax,title in zip(axs,['Tracking integrator controls','Fluid and saved-data sensitivity']):ax.axhline(threshold,ls=':',color='.3',label='5% initial spacing screen');ax.set(title=title,xlabel='Model time',ylabel='Maximum marker position difference');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Passing a tracking-step check does not establish fluid-grid agreement');fig.savefig(ROOT/'sensitivity.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for g,ax in enumerate(axs):
        for n,style,label in [('reference','-','64 grid'),('grid128_refined','--','128 grid')]:
            r=runs[n];ax.semilogy(r['times'],[v[g]['seed_amplification_rms'] for v in r['metrics']],style,label=label)
        ax.set(title=titles[g],xlabel='Model time',ylabel='RMS separation / initial displacement');ax.legend();ax.grid(alpha=.2)
    fig.suptitle('Small position changes grow differently across grids | finite-time sensitivity');fig.savefig(ROOT/'position-sensitivity.png',dpi=150);plt.close(fig)
    fig=plt.figure(figsize=(11,5),layout='constrained');z=raw['reference'];p=z['positions']
    for g in (0,1):
        ax=fig.add_subplot(1,2,g+1,projection='3d')
        for i in range(g*125,(g+1)*125):ax.plot(*p[:,i].T,color=colors[g],alpha=.16,lw=.6)
        ax.scatter(*p[0,g*125:(g+1)*125].T,s=9,c=colors[g],label='Start');ax.scatter(*p[-1,g*125:(g+1)*125].T,s=10,c='black',label='End')
        ax.set(title=titles[g],xlabel='x (unwrapped)',ylabel='y (unwrapped)',zlabel='z (unwrapped)');ax.legend(fontsize=8)
    fig.suptitle('125 paths per neighborhood in the 64-grid reconstruction | diagnostic trajectories');fig.savefig(ROOT/'trajectories.png',dpi=150);plt.close(fig)
    report=['# Individual fluid markers: stress test','','**The markers develop different motion. Their detailed paths fail the grid and saved-time sensitivity checks.**','','The user selected individual fluid markers following the approved Navier–Stokes flow. This completed test reconstructs passive tracks from existing velocity fields through model time 0.30. It does not add a molecular model, particle forces, noise or a new fluid run. The ongoing 128-grid run contributes only already-saved observations through 0.30; its full run is not claimed complete.','','## What individuality means here','','Each label follows `dX_i/dt = u_h(X_i,t)`. Different positions sample different velocities. Within a cloud, `c_i = u_h(X_i,t) - mean_j u_h(X_j,t)` measures the difference from that cloud’s sampled mean. This c is a diagnostic of variations in the resolved/interpolated flow, not molecular thermal motion. The cloud mean is an equally weighted marker average, not a mass-weighted average of a filled fluid volume. Five exact copies of initial positions follow identical paths: a different label alone does not change the dynamics.','','## Stress and controls','','- 125 labels fill a 5x5x5 central neighborhood, and 125 fill a neighborhood of the original 64-grid initial global spin peak at (-2.90625,-2.90625,-2.90625), crossing the periodic join. Initial nearest lattice spacing is 0.09375; cloud side length is 0.375.','- 600 fixed links connect adjacent lattice sites. A separate nearest-six-neighbor diagnostic measures changes in proximity, with stable label-index tie breaking.','- Each label gets a companion initially displaced by exactly 0.0009375 (1% of spacing), using random-direction seed 104. This changes marker positions only.','- Eight reconstructions test smaller RK4 tracking steps, a halved fluid timestep, finer fluid grid and coarser saved-time spacing. The additional eighth check investigates the failed initial 128-grid tracking-step screen.','- The fluid remains the supplied original start, viscosity 0.001, periodic cube side 6, Heun evolution, preserved per-grid mean and zero external force.','','## Measured group behavior','','All positions, speeds and times use the original simulation model units.','','![Individual motion](individual-motion.png)','','| Grid | Neighborhood | Final median original-neighbor spacing | Final RMS velocity difference from cloud mean | Final initial-neighbor retention |','| --- | --- | ---: | ---: | ---: |']
    for n,grid in [('reference',64),('grid128_refined',128)]:
        for g in (0,1):
            m=runs[n]['metrics'][-1][g];report.append(f"| {grid} | {titles[g]} | {m['neighbor_distance_quantiles'][2]:.5f} | {m['velocity_deviation_rms']:.4f} | {100*m['initial_neighbor_retention']:.1f}% |")
    report+=['','Every link starts with distance 0.09375. These values describe the sampled trajectories; differences between the rows are part of the result. Approach or separation follows the given flow and does not demonstrate a new attraction law. Crowding of labels is not fluid-density increase.','','![Trajectories](trajectories.png)','','Coordinates above are continuously unwrapped across the periodic box. They are three-dimensional physical paths, not an autonomous phase portrait. Projected crossings and sparse sampled near approaches do not demonstrate molecule collisions or disprove trajectory uniqueness.','','## Did the paths survive the controls?','','The declared exploratory screen is maximum corresponding-label displacement no greater than 0.0046875 (5% of initial spacing). It is a descriptive tolerance, not a rigorous physical error bound.','','| Comparison | Maximum through 0.10 | Maximum through 0.30 | Screen through 0.30 |','| --- | ---: | ---: | --- |']
    for label,c in comparisons.items():report.append(f"| {label} | {c['maximum_through_0_1']:.7f} | {c['maximum_through_0_3']:.7f} | {'PASS' if c['screen_through_0_3'] else 'FAIL'} |")
    grid=comparisons['64 to 128 grid'];ref=comparisons['128 tracking step: fine to refined']
    report+=['','![Numerical controls](sensitivity.png)','',f"The extra 128-grid half-step changes the maximum tracked position by {ref['maximum_through_0_3']:.6f}. The 64/128 grid comparison changes it by {grid['maximum_through_0_3']:.6f}. Subtracting the difference in each grid’s preserved uniform mean drift still leaves a maximum discrepancy of {grid['mean_drift_adjusted_maximum']:.6f}. This frame adjustment is diagnostic only; it does not equalize the starting fields or their evolution.",'','The grids start with different sampled fields, mean velocities and energies. Their discrepancy therefore combines starting-field and spatial-discretization differences; it is not a clean same-smooth-start convergence result. Temporal coarsening from 0.01 to 0.02 is also a sensitivity test, not an upper bound on the error at0.01. The raw periodic-join mismatch and unresolved vorticity peaks remain.','','## Small starting-position changes','','![Position sensitivity](position-sensitivity.png)','','| Grid | Neighborhood | Final RMS position amplification | Largest final individual amplification |','| --- | --- | ---: | ---: |']
    for n,grid in [('reference',64),('grid128_refined',128)]:
        for g in (0,1):
            m=runs[n]['metrics'][-1][g];report.append(f"| {grid} | {titles[g]} | {m['seed_amplification_rms']:.3f} | {m['seed_amplification_max']:.3f} |")
    report+=['','Amplification is shortest-periodic companion separation divided by its initial displacement. These are finite-time marker sensitivity ratios, not asymptotic Lyapunov exponents. A single perturbation direction per label and these two grids do not establish a universal rate, chaos classification, periodic breathing or physical molecular behavior.','','## Verification and reproduction','','All 93 distinct source fields were finite-checked and hashed. Initial velocity interpolation was cross-checked against SciPy; a known rigid-rotation trajectory checks RK4 convergence. Every stored trajectory is finite, preserves identical-clone agreement and has its own hash. Final sampled velocities and reported group diagnostics are checked by the separate verification script. Numerical fluid sources remain unchanged.','','Rebuild this report and four charts with `python report.py` beside the supplied JSON and NPZ files. `run.py` and `refine.py` require the original workspace’s saved fluid fields. Reproduction does not require editing or rerunning the fluid solver.','','[Changing separation roles](separation-roles.md) · [Fractal tests](fractal/README.md) · [Protocol](protocol.json) · [Results](results.json) · [Comparison details](summary.json) · [Verification](verification.json) · [Earlier tracked pairs](../central-response/tracked-patterns/README.md) · [Main study](../README.md)','']
    (ROOT/'README.md').write_text('\n'.join(report),encoding='utf-8');print(json.dumps({'comparisons':{k:{j:v[j] for j in ['maximum_through_0_3','screen_through_0_3']} for k,v in comparisons.items()},'fields':len(data['fields'])},indent=2))
if __name__=='__main__':main()
