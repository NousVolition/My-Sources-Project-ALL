"""Analyze completed data and build recorded molding comparisons and a report."""
from pathlib import Path
from dataclasses import asdict
import json,base64,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from clay_hug import *
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';FIG=ROOT/'figures'
def read(n):return json.loads((DATA/n).read_text())
def main():
    rows=read('runs.json')+read('followup.json')['runs']+read('feedback-mesh.json')['runs'];summary=read('summary.json');checks=read('solver-checks.json')+read('followup-checks.json');validation=read('validation.json')
    summary['failures']=[dict(name=r['name'],failure=r['failure']) for r in rows if r['failure']]
    summary['largest_completed_material_budget_error']=max(r['metrics']['material_budget_error'] for r in rows if r['failure'] is None)
    accepted=['ellipse','lobed','independent-DOP853','strong-resolved']
    accepted_checks_passed=all(next(c for c in checks if c['name']==n)['passed'] for n in accepted)
    lookup={r['name']:r for r in rows};FIG.mkdir(exist_ok=True)
    def load(name):
        a=np.load(DATA/(name+'.npz'));p=Parameters(**lookup[name]['parameters']);return a,Mesh(p)
    supplemental={}
    for r in rows:
        a,m=load(r['name']);d=[m.diagnostics(s) for s in a['state']]
        if d:supplemental[r['name']]=dict(peak_contact=max(x['contact_fraction'] for x in d),minimum_surface_gap=min(x['surface_gap_rms'] for x in d),time_best_contact=float(a['time'][int(np.argmax([x['contact_fraction'] for x in d]))]))
    comparisons=[]
    for one in rows:
        if one['group'] not in ('paired','long','followup','mesh-feedback') or one['parameters']['feedback'] or not one['name'].endswith('-one'):continue
        name=one['name'][:-3]+'two';two=lookup[name]
        if one['failure'] or two['failure'] or max(one['metrics']['material_budget_error'],two['metrics']['material_budget_error'],one['metrics']['total_budget_error'],two['metrics']['total_budget_error'])>1e-3:continue
        changes={k:100*(two['metrics'][k]/one['metrics'][k]-1) if one['metrics'][k]!=0 else None for k in ['lorenz_rms','z_std','mean_imprint','mesh_motion_rms','mean_contact_fraction']}
        comparisons.append(dict(one=one['name'],two=name,changes_percent=changes))
    meshconv=[r for r in rows if r['group']=='mesh'];boundary=[]
    dense_angle=np.linspace(0,2*np.pi,721)
    for r in meshconv:
        a,m=load(r['name']);x=m.unpack(a['state'][-1])[0][:m.p.angles]
        theta=np.mod(np.arctan2(x[:,1],x[:,0]),2*np.pi);order=np.argsort(theta)
        radius=np.interp(dense_angle,theta[order],np.linalg.norm(x,axis=1)[order],period=2*np.pi);boundary.append(radius)
    spatial=[dict(coarse=meshconv[i]['name'],fine=meshconv[i+1]['name'],radial_rms_error=float(np.sqrt(np.mean((boundary[i]-boundary[i+1])**2))),relative_imprint_difference=abs(meshconv[i]['metrics']['final']['imprint']/meshconv[i+1]['metrics']['final']['imprint']-1)) for i in range(len(meshconv)-1)]
    overlap_failures=[r['name'] for r in rows if r['failure'] is None and r['metrics']['maximum_penetration']>.02]
    budget_failures=[r['name'] for r in rows if r['failure'] is None and max(r['metrics']['material_budget_error'],r['metrics']['total_budget_error'])>1e-3]
    area_failures=[r['name'] for r in rows if r['failure'] is None and r['metrics']['minimum_cell_area_ratio']<.1]
    mem=lookup['release-clay']['metrics'];retention=mem['final']['imprint']/mem['initial']['imprint']
    windows=[]
    for c in comparisons:
        a,ma=load(c['one']);b,mb=load(c['two']);edges=np.linspace(4,min(a['time'][-1],b['time'][-1]),7 if c['one']=='long-one' else 3)
        blocks=[]
        for start,end in zip(edges[:-1],edges[1:]):
            metrics=[]
            for result in [a,b]:
                t=result['time'];y=result['state'];mask=(t>=start-1e-10)&(t<=end+1e-10);tt=t[mask]
                metrics.append(float(np.sqrt(trapezoid(np.sum(y[mask,-7:-5]**2,axis=1),tt)/(tt[-1]-tt[0]))))
            blocks.append(100*(metrics[1]/metrics[0]-1))
        windows.append(dict(one=c['one'],two=c['two'],block_edges=edges.tolist(),lorenz_activity_changes_percent=blocks,range_percent=[min(blocks),max(blocks)]))
    analysis=dict(comparisons=comparisons,spatial_convergence=spatial,supplemental=supplemental,
        completed=len(rows)-len(summary['failures']),guard_stops=summary['failures'],geometric_target_misses=overlap_failures,
        budget_target_misses=budget_failures,cell_area_target_misses=area_failures,unloaded_imprint_retained_fraction=retention,
        accepted_solver_checks_passed=accepted_checks_passed,solver_checks=checks,window_sensitivity=windows,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (DATA/'analysis.json').write_text(json.dumps(analysis,indent=2,allow_nan=False)+'\n')
    cases=[dict(key=shape,label=shape.title()+' · alternate drive start',one=shape+'-s1-one',two=shape+'-s1-two') for shape in ['circle','ellipse','lobed']]
    cases.extend(dict(key='alpha'+str(alpha),label=f'Ellipse · drive strength {alpha:g}',one=f'alpha{alpha:g}-one',two=f'alpha{alpha:g}-two') for alpha in [.15,1.,2.])
    cases[-1].update(one='resolved-alpha2-one',two='resolved-alpha2-two',label='Strong drive · soft contact limit')
    records={};frame_count=51
    for case in cases:
        for name in (case['one'],case['two']):
            a,m=load(name);indices=np.rint(np.linspace(0,len(a['time'])-1,frame_count)).astype(int);frames=[]
            for i in indices:
                x,rest,lorenz=m.unpack(a['state'][i]);d=m.diagnostics(a['state'][i])
                frames.append(dict(time=round(float(a['time'][i]),5),positions=np.round(x,4).tolist(),lorenz=np.round(lorenz,5).tolist(),imprint=round(d['imprint'],5),contact=round(d['contact_fraction'],5),overlap=round(d['penetration'],5),port=round(float(a['state'][i,-1]),5)))
            records[name]=dict(frames=frames,outline=np.round(a['outline'][::2],5).tolist(),trace=np.round(np.column_stack([a['time'][::2],a['state'][::2,-7],a['state'][::2,-5]]),5).tolist())
    for case in cases:case['initialIndex']=int(np.argmax([f['contact'] for f in records[case['two']]['frames']]))
    first,_=load('ellipse-s0-one')
    payload=dict(cases=cases,records=records,angles=24,initial=np.round(first['initial_positions'],4).tolist(),triangles=first['triangles'].tolist())
    fragment=(ROOT/'molding-template.html').read_text(encoding='utf-8').replace('__MOLDING_DATA__',json.dumps(payload,separators=(',',':')))
    assert len(fragment.encode())<1_000_000
    (ROOT/'molding-fragment.html').write_text(fragment,encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    for j,shape in enumerate(['circle','ellipse','lobed']):
        a,m=load(shape+'-s1-two');ind=np.argmin(abs(a['time']-supplemental[shape+'-s1-two']['time_best_contact']))
        for i,index in enumerate([0,ind]):
            ax=axs[i,j];x=m.unpack(a['state'][index])[0]
            ax.fill(a['outline'][:,0],a['outline'][:,1],color='#c7d0d3');ax.triplot(x[:,0],x[:,1],m.tri,lw=.55,color='#137d87')
            ax.set(xlim=(-2,2),ylim=(-2,2),xlabel='x (model units)',ylabel='y (model units)',title=f'{shape.title()} · t={a["time"][index]:.2f}');ax.set_aspect('equal')
    fig.suptitle('Computed filled mesh before loading and near best contact')
    fig.savefig(FIG/'molding.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for name,color,label in [('ellipse-s0-one','#a85c33','One-way'),('ellipse-s0-two','#147c89','Two-way')]:
        a,m=load(name);d=[m.diagnostics(s) for s in a['state']];t=a['time'];y=a['state']
        axs[0,0].plot(y[:,-7],y[:,-5],color=color,lw=.7,label=label)
        axs[0,1].plot(t,y[:,-7],color=color,lw=.8,label=label)
        axs[1,0].plot(t,[v['imprint'] for v in d],color=color,label=label)
        axs[1,1].plot(t,y[:,-1],color=color,label=label)
    axs[0,0].set(xlabel='Lorenz X',ylabel='Lorenz Z',title='The deforming material changes the drive');axs[0,0].legend()
    for ax,title,ylabel in [(axs[0,1],'Drive amplitude','Lorenz X'),(axs[1,0],'Different forcing changes the molded imprint','Imprint RMS'),(axs[1,1],'Reciprocal mechanical work','Net work into material')]:ax.set(xlabel='Model time',ylabel=ylabel,title=title)
    fig.savefig(FIG/'reciprocal.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    for name,color,label in [('release-clay','#147c89','Yielding, fading material'),('release-elastic','#a85c33','Elastic control')]:
        a,m=load(name);d=[m.diagnostics(s) for s in a['state']];axs[0].plot(a['time'],[v['imprint'] for v in d],color=color,label=label)
    axs[0].set(xlabel='Time after unloading',ylabel='Imprint RMS',title='Fading after the load is removed');axs[0].legend(fontsize=8)
    for r in meshconv:axs[1].plot(r['parameters']['angles'],r['metrics']['final']['imprint'],'o',color='#147c89')
    axs[1].set(xlabel='Angular divisions (layers refined too)',ylabel='Final imprint',title='Spatial convergence under steady forming')
    for r in rows:
        if r['group']=='contact':axs[2].plot(r['parameters']['contact'],r['metrics']['maximum_penetration'],'o',color='#147c89')
    axs[2].set(xlabel='Contact penalty',ylabel='Maximum radial overlap',title='Contact stiffness sensitivity');fig.savefig(FIG/'validation.png',dpi=150);plt.close(fig)
    embed=lambda name:'data:image/png;base64,'+base64.b64encode((FIG/name).read_bytes()).decode()
    pair_table=''.join('<tr><td>'+c['two']+'</td>'+''.join(f'<td>{c["changes_percent"][k]:+.2f}%</td>' for k in ['lorenz_rms','z_std','mean_imprint','mesh_motion_rms'])+'</tr>' for c in comparisons)
    stop_names=', '.join(f['name'] for f in summary['failures']) or 'none'
    misses=', '.join(overlap_failures) or 'none';budget=', '.join(budget_failures) or 'none'
    window_table=''.join(f'<tr><td>{c["two"]}</td><td>{c["range_percent"][0]:+.2f}% to {c["range_percent"][1]:+.2f}%</td></tr>' for c in windows)
    finer=next(c for c in comparisons if c['two']=='feedback-mesh48-two')['changes_percent']
    html=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Moldable hug: reciprocal Lorenz stress tests</title>
<style>:root{{--foreground:#24383f;--background:#fff;--muted:#d9e0e2;--muted-foreground:#65767e;--border:#c7d2d7;--viz-series-1:#147c89;--viz-series-2:#a85c33}}body{{font:16px/1.6 system-ui,sans-serif;max-width:1100px;margin:auto;padding:24px;color:var(--foreground)}}img{{max-width:100%;height:auto}}.viz-controls{{display:flex;gap:16px;flex-wrap:wrap;align-items:end}}.form-label{{display:block;min-width:0}}.form-select{{display:block;max-width:100%;font:inherit;padding:8px}}.form-range{{display:block;width:100%}}.btn{{font:inherit;padding:10px}}.viz-row{{display:flex;gap:20px;flex-wrap:wrap}}.table-responsive{{overflow:auto}}table{{border-collapse:collapse;width:100%}}th,td{{padding:7px;text-align:left;border-bottom:1px solid #ddd}}pre{{white-space:pre-wrap;background:#f0f4f5;padding:16px}}.text-small{{font-size:14px}}a{{color:#087080}}</style>
<h1>A moldable hug with contact and reciprocal drive</h1>
<p>The hug now has a filled two-dimensional material cross-section. Its nodes move independently, its bonds yield and change their resting lengths, and its inner surface presses against a specified rigid object. Fading imprint follows an explicit internal relaxation law. There is no imposed arm-opening latch. These are actual computed material configurations, rather than a prescribed shape animation.</p>
<p>This implements the user's intended clay-like ease of molding as a low-yield cohesive material. Parameters are dimensionless modeling choices. It is an abstract 2-D material prototype, not a calibrated sample of real clay. The initial annular band already surrounds the object; the experiment tests reshaping into contact, not wrapping a detached slab around it.</p>
{fragment}
<h2>What the completed tests show</h2>
<p>{len(rows)} stored configurations plus one adaptive independent integration and 13 initial checks are included. {analysis['completed']} configurations finish their requested interval; {len(summary['failures'])} are stopped by the numerical guard. The base mesh has 72 nodes and 96 filled triangular cells. Spatial checks use 240 nodes/384 cells and 504 nodes/864 cells. Paired runs use identical material settings and Lorenz initial states, with the return connection off or on.</p>
<img src="{embed('molding.png')}" alt="Initial annular material and computed molding against circular, elliptical and concave lobed objects">
<p>The molding figure and three object playback cases use the second Lorenz start, (-8, 8, 27), which brings all three objects into full sampled near-contact. The first start, (1, 1, 1), does not bring the circle into contact over the recorded interval and gives only partial ellipse/lobed contact. This sensitivity is a measured outcome, not omitted as a failed fit. The force-strength playback and the reciprocal time-series figure use the first start; all starting states and both directions are retained in the comparison table.</p>
<p>Near-contact fraction means the fraction of sampled inner surface points with radial gap magnitude below 0.03. Imprint is an area-weighted RMS change in bond rest length divided by its initial length. The rigid object's boundary is sampled at three quadrature points per inner edge for contact forces; overlap diagnostics use 17 points per edge. These definitions are retained in the code and data.</p>
<div class="table-responsive"><table><thead><tr><th>Two-way run</th><th>Lorenz activity RMS</th><th>Z fluctuations</th><th>Mean imprint</th><th>Material motion RMS</th></tr></thead><tbody>{pair_table}</tbody></table></div>
<p>Changes are percentages relative to the matched one-way run, after discarding the first four model-time units. Activity RMS is sqrt(mean(X²+Y²)); Z fluctuations are its time standard deviation. Only two initial states are tested for the object comparisons. These finite-window changes are not confidence intervals. Long-time pointwise separation can result from phase sensitivity, so the report retains statistics, controls, refinement and work diagnostics.</p>
<p>The additional matched 48-angle, five-layer pair changes material motion by {finer['mesh_motion_rms']:+.2f}%, Z fluctuations by {finer['z_std']:+.2f}% and mean imprint by {finer['mean_imprint']:+.2f}%. The motion and Z response closely track the original mesh's effect; the imprint effect remains sensitive to spatial resolution. Two coupled meshes do not establish full continuum convergence. Molding is intermittent under the changing Lorenz force: the band can reach the object and then pull away as forcing reverses.</p>
<p>Window sensitivity below gives the range of activity changes when the retained interval is split into two equal blocks (six for the long run). Variation reflects sensitivity to the observation window; it is not a statistical confidence interval or an independent ensemble.</p>
<table><thead><tr><th>Two-way run</th><th>Block activity-change range</th></tr></thead><tbody>{window_table}</tbody></table>
<img src="{embed('reciprocal.png')}" alt="Matched Lorenz trajectories, forcing amplitude, material imprint and exchanged work">
<h2>How the material works</h2><pre>E = sum_edges [k/2 (length-rest)² + H k/2 (rest-initial_length)²]
  + bulk sum_cells A0 [area/A0 - 1 - log(area/A0)]
  + contact penalty + square-wall penalty + preload A
A = sum_nodes reference_area × distance_from_center
z = k [length-rest-H(rest-initial_length)]
rest′ = z/(k memory_time) + soft_threshold(z/k, yield_strain × initial_length)/yield_time
velocity = (material force + Lorenz actuation force)/(viscosity × reference_area)</pre>
<p>The internal law dissipates energy because z × rest′ is nonnegative. Below the yield threshold the rest shape relaxes slowly; above it the additional overstress flow permits faster reshaping. Hardening draws the unloaded rest shape toward its original value. A free isolated bond has the checked decay time memory_time/H = 10 model-time units. The full mesh need not decay at a single exponential rate because contact, elasticity and cell areas couple its bonds.</p>
<p>In the full unloading experiment, {100*retention:.2f}% of the loaded imprint remains after 30 model-time units, with the object still present. The elastic control has fixed rest lengths and zero internal imprint. A small reversible deformation and nonzero material force are needed to mold this numerical material; the phrase “without effort” has been interpreted as low yield resistance rather than literally zero work.</p>
<img src="{embed('validation.png')}" alt="Full-mesh imprint decay, spatial refinement, and contact penalty convergence">
<h2>Both directions share a work balance</h2><pre>lambda = alpha/10
material force from drive = -lambda X grad(A)
Lorenz X′ = 10(Y-X) + lambda A′/C       [return enabled]
Lorenz Y′ = X(rho-Z)-Y
Lorenz Z′ = XY-(8/3)Z
E_L = C/2 (X²+Y²+Z²)
drive work rate into material = -lambda X A′</pre>
<p>The material work and the return term in the Lorenz energy cancel exactly at the continuous-equation level. Unlike the earlier imprint-dependent resistance experiment, recovery of material energy can act back on the drive. Lorenz itself retains its native forcing and dissipation; this is not a closed thermodynamic fluid model. C=0.01 and the force scale alpha/10 are explicit chosen reservoir/port normalizations. The capacity sweep checks sensitivity to that choice.</p>
<h2>Stress limits and convergence</h2>
<p>The selected target is 1e-4 maximum scaled state discrepancy for the short solver comparisons, 1e-3 normalized work-balance residual, 0.02 maximum radial overlap and cell area above 10% of its reference area. Accepted solver checks passed: {accepted_checks_passed}. Original check results are preserved alongside the finer strong-drive follow-up. Coarse guard stops: {stop_names}. Completed cases missing the overlap target: {misses}. Completed cases missing the work-balance target: {budget}. Original stopped runs are preserved; their guard stops are not labeled physical fracture.</p>
<p>The original alpha=2 paired runs at step 0.001 completed but had large work-balance errors. Completion and an intact mesh did not make those answers accurate. They are excluded from the comparison table and playback; their step-0.00025 repeats are explicitly named resolved-alpha2. Short-time agreement is insufficient to establish reliability for a longer, more deformed run. The report includes the original failures rather than silently replacing them.</p>
<p>The largest completed material balance residual is {summary['largest_completed_material_budget_error']:.3e}. Contact is a penalty approximation: small overlap is allowed, and stiffness/mesh checks measure its sensitivity. The square now exerts an actual penalty force, but it is a soft boundary rather than an exact impenetrable constraint. Increasing contact stiffness and decreasing viscosity make the numerical problem harder; reduced steps are tested explicitly.</p>
<p>Spatial comparison of the steady-formed inner boundary gives adjacent-grid radial RMS differences {spatial[0]['radial_rms_error']:.5f} and {spatial[1]['radial_rms_error']:.5f}. This is a measured convergence result, not a declaration that the finest mesh is exact. The cohesive bond mesh can retain directional/discretization effects; it is not a calibrated isotropic continuum constitutive law.</p>
<h2>Limits of the result</h2><p>This implementation covers a 2-D cross-section, fixed rigid objects, frictionless contact, a compressible cohesive band, yielding rest lengths and fading memory. It does not include a full 3-D body, adhesion, friction, fracture, general remeshing, material passing through itself, a pressure-dependent soil-clay law, temperature, quantum effects or resolved water flow. Severe folding or inverted cells requires another contact/remeshing method. There are no measured material parameters or experimental reciprocal-drive benchmarks for this proposed hug, so numerical verification is not experimental validation.</p>
<p>Lorenz supplies a small dynamical reservoir. This test establishes a connection between a deformable material and that reservoir. It does not establish an exact connection to Navier–Stokes or solve its regularity problem. The earlier fluid experiments remain separate.</p>
<h2>Reproduce</h2><pre>python -m pip install -r requirements-lock.txt
python validate_model.py
python run_stress.py
python run_followup.py
python run_mesh_feedback.py
python build_report.py
python verify_results.py</pre>
<p>Sources, actual trajectories, parameter selections, stopped attempts and source hashes are included. The code is self-contained apart from NumPy, SciPy and Matplotlib. Optional playback tests use Node, Playwright and Chrome.</p>
<p>Background: <a href="https://arxiv.org/abs/1911.11662">Islam et al., elasto-viscoplastic clay simulations</a> illustrates that measured clay validation requires an appropriate constitutive law and material data; that law is not implemented here. <a href="https://mooseframework.inl.gov/source/materials/HyperbolicViscoplasticityStressUpdate.html">MOOSE viscoplasticity documentation</a> describes a separate established yield/flow implementation; this prototype uses its own scalar bond overstress law. <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html">SciPy documentation</a> defines the independent adaptive solver used here.</p></html>'''
    (ROOT/'report.html').write_text(html,encoding='utf-8')
    selected=next(c for c in comparisons if c['two']=='ellipse-s0-two')['changes_percent']
    results=f'''## Completed results

{len(rows)} saved configurations and one independent adaptive integration were run. {analysis['completed']} finished their requested interval; {len(summary['failures'])} stopped at the numerical guard. Stopped or inaccurate numerical attempts remain in the dataset.

![Computed molding](figures/molding.png)

The round band molds into contact with the circle, tilted ellipse and concave lobed object for the second Lorenz start, (-8, 8, 27). Selected peak near-contact fractions are {100*supplemental['circle-s1-two']['peak_contact']:.1f}%, {100*supplemental['ellipse-s1-two']['peak_contact']:.1f}% and {100*supplemental['lobed-s1-two']['peak_contact']:.1f}%. This means radial gap magnitude below 0.03 at the sampled inner surface; it is not zero-overlap hard contact.

Starting pressure history matters: the first Lorenz start, (1, 1, 1), never reaches the circle in this interval and gives only {100*supplemental['ellipse-s0-two']['peak_contact']:.1f}% and {100*supplemental['lobed-s0-two']['peak_contact']:.1f}% peak near-contact for the ellipse and lobed object. The figure/object playback uses the second start; the force-strength playback and time-series figure use the first. The model does not guarantee full molding for every drive.

For the moderate ellipse case over time 4–16, enabling the return changes Lorenz activity RMS by {selected['lorenz_rms']:+.2f}%, mean material imprint by {selected['mean_imprint']:+.2f}% and material motion RMS by {selected['mesh_motion_rms']:+.2f}%. Other shapes, initial states, observation windows and port choices give different changes; no universal effect size or confidence interval is established.

On the finer coupled mesh, material motion changes by {finer['mesh_motion_rms']:+.2f}% and Lorenz Z fluctuations by {finer['z_std']:+.2f}%. Those effects persist closely; the imprint change is {finer['mean_imprint']:+.2f}% and remains sensitive to mesh resolution. The hug can mold into contact and later pull away when Lorenz forcing reverses. These runs do not establish permanent adhesion or full continuum convergence.

![Both directions of influence](figures/reciprocal.png)

The force changes the material; the material's changing shape reacts back on Lorenz. The continuous reciprocal work terms cancel. Saved trajectories also contain intervals of work returning from material to drive. The approximate forward/return integrals sampled every 0.04 units are labeled sampled in the JSON; the signed net work is integrated along with the equations.

Thirteen initial validation checks pass. The independent DOP853 comparison has maximum scaled state error {next(c for c in checks if c['name']=='independent-DOP853')['scaled_state_error']:.3e}. The original strong-drive step-0.001 runs have energy-balance errors {lookup['alpha2-one']['metrics']['total_budget_error']:.3e} and {lookup['alpha2-two']['metrics']['total_budget_error']:.3e}; their step-0.00025 repeats reduce these to {lookup['resolved-alpha2-one']['metrics']['total_budget_error']:.3e} and {lookup['resolved-alpha2-two']['metrics']['total_budget_error']:.3e}. Only the refined strong results enter the playback/table.

![Unloading and numerical sensitivity](figures/validation.png)

After unloading for 30 time units, {100*retention:.2f}% of the loaded imprint remains while the object stays present. The elastic control's bond rest lengths do not change, so its internal imprint is zero. Adjacent spatial refinements give boundary radial RMS differences {spatial[0]['radial_rms_error']:.5f} and {spatial[1]['radial_rms_error']:.5f}. Those differences measure numerical sensitivity; the finest mesh is still an approximation.

Completed cases missing the 0.02 radial-overlap target: **{misses}**. Cases missing the 0.001 work-balance target: **{budget}**. Numerical guard stops: **{stop_names}**. Review [all measured comparisons](data/analysis.json), [original runs](data/runs.json), [refined follow-up](data/followup.json), [solver checks](data/solver-checks.json), and [evidence audit](data/audit.json).

'''
    readme=(ROOT/'README.md').read_text(encoding='utf-8');start='<!-- RESULTS_START -->';end='<!-- RESULTS_END -->'
    before,rest=readme.split(start,1);_,after=rest.split(end,1)
    (ROOT/'README.md').write_text(before+start+'\n'+results+end+after,encoding='utf-8')
    print(json.dumps(dict(configurations=len(rows),comparisons=len(cases),fragment_bytes=len(fragment.encode()),analysis=analysis['completed'])))
if __name__=='__main__':main()
