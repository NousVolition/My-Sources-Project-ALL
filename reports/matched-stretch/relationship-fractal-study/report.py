"""Build the local report and scientific figures from saved analysis."""
from pathlib import Path
import json, html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
def read(n):return json.loads((ROOT/n).read_text(encoding='utf-8'))
def main():
    A=read('analysis.json');M=A['metrics'];C=read('results.json')['cases'];S=A['sensitivity']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fafaf7','axes.facecolor':'#fafaf7','savefig.facecolor':'#fafaf7'})
    names=['menger','koch','mixed','circle','ball'];titles=['Menger sample','Koch boundary','Menger + Koch','Circle control','Ball control']
    colors=['#126a76','#d18b22','#8756ad','#466cc2','#62724b']
    profiles=['fading','mutual_amplification','opposing_responses'];ptitles=['Fading','Mutual amplification','Opposing responses'];initials=['balanced','positive','negative']
    def save(fig,n):
        fig.savefig(ROOT/(n+'.png'),dpi=165,bbox_inches='tight');fig.savefig(ROOT/(n+'.svg'),bbox_inches='tight');plt.close(fig)
    con=np.load(ROOT/'construction.npz');fig=plt.figure(figsize=(15,3.9),layout='constrained')
    for k,(name,title,color) in enumerate(zip(names,titles,colors)):
        ax=fig.add_subplot(1,5,k+1,projection='3d');x=con[name];ax.scatter(*x.T,s=9,c=color,alpha=.8,depthshade=False)
        ax.set(title=title,xlabel='x',ylabel='y',zlabel='z',xlim=(-.3,.3),ylim=(-.3,.3),zlim=(-.3,.3));ax.set_box_aspect((1,1,1));ax.set_xticks([-.2,0,.2]);ax.set_yticks([-.2,0,.2]);ax.set_zticks([-.2,0,.2])
    fig.suptitle('Five starting arrangements • 192 markers each • RMS radius 0.15',fontsize=16);save(fig,'01-geometries')
    fig,axes=plt.subplots(1,2,figsize=(12,7.2),layout='constrained');arrays=[]
    for run in ['step64','grid128']:
        pairs=M[run]['pairs'];arrays.append(np.array([[next(p['radius_change_percent'] for p in pairs if p['geometry']==n and p['profile']==pr and p['initial_state']==st) for n in names] for pr in profiles for st in initials]))
    bound=max(abs(a).max() for a in arrays)
    for ax,arr,title in zip(axes,arrays,['64-grid saved flow','128-grid saved flow']):
        im=ax.imshow(arr,cmap='RdBu_r',vmin=-bound,vmax=bound,aspect='auto')
        ax.set_xticks(range(5),['Menger','Koch','Mixed','Circle','Ball']);ax.set_yticks(range(9),[p+' / '+s for p in ptitles for s in initials]);ax.set_title(title)
        for i in range(9):
            for j in range(5):ax.text(j,i,f'{arr[i,j]:+.1f}',ha='center',va='center',color='white' if abs(arr[i,j])>.56*bound else '#15242e',fontsize=10)
        for y in [2.5,5.5]:ax.axhline(y,color='white',lw=2)
    fig.colorbar(im,ax=axes,label='Cloud RMS radius change versus passive markers (%)',shrink=.85)
    fig.suptitle('Relationships can contract or expand the cloud in this chosen model\nAt t = 0.10; negative = smaller than passive, positive = larger',fontsize=15);save(fig,'02-radius-effects')
    d=np.load(ROOT/'step64.npz');fig,axes=plt.subplots(1,3,figsize=(12,4.6),layout='constrained')
    ids=[next(i for i,c in enumerate(C) if c['geometry']=='mixed' and c['profile']=='mutual_amplification' and c['initial_state']==s and c['mu']==mu) for s,mu in [('positive',0),('positive',20),('negative',20)]]
    allx=d['positions'][:,ids];lo=allx[:,:,:,0].min();hi=allx[:,:,:,0].max();bottom=allx[:,:,:,1].min();top=allx[:,:,:,1].max()
    for ax,idx,title in zip(axes,ids,['Passive markers','Positive states + added attraction','Negative states + added repulsion']):
        x=d['positions'][:,idx]
        for j in range(0,192,6):ax.plot(x[:,j,0],x[:,j,1],color='#9b9b98',lw=.7,alpha=.5)
        ax.scatter(x[0,:,0],x[0,:,1],s=8,c='#aeb5b7',alpha=.5,label='Start');ax.scatter(x[-1,:96,0],x[-1,:96,1],s=13,c=colors[0],label='Menger labels');ax.scatter(x[-1,96:,0],x[-1,96:,1],s=13,c=colors[1],label='Koch labels')
        ax.set(title=title,xlabel='x',ylabel='y',xlim=(lo-.02,hi+.02),ylim=(bottom-.02,top+.02));ax.set_aspect('equal',adjustable='box')
    axes[0].legend(fontsize=8,loc='upper left');fig.suptitle('The mixed cloud under three motion conditions\nMutual-amplification profile, 64-grid flow; x–y projections can overlap',fontsize=15);save(fig,'03-mixed-motion')
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    for ax,pr,title in zip(axes,profiles,ptitles):
        for name,color in zip(names,colors):
            for mu,style in [(0,'--'),(20,'-')]:
                c=next(c for c in M['step64']['cases'] if c['geometry']==name and c['profile']==pr and c['initial_state']=='positive' and c['mu']==mu)
                ax.plot(M['step64']['times'],np.array(c['state_rms_curve'])/c['state_rms_initial'],style,color=color,label=name if mu else None,lw=1.8)
        ax.set(title=title,xlabel='Time',ylabel='State RMS / initial RMS');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('Relationship response depends on the chosen coefficients\nPositive initial states • solid: added motion • dashed: passive',fontsize=15);save(fig,'04-state-response')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for label,color in [('step','#126a76'),('snapshot','#d18b22'),('grid','#8756ad')]:
        for ax,key in zip(axes,['max_position_difference','max_relative_state_rms_difference']):
            val=np.sort([r[key] for r in S[label]['cases']]);ax.semilogy(range(1,91),np.maximum(val,1e-17),label=label,color=color,lw=2)
    axes[0].axhline(.001,color='#777',ls=':',label='Step-screen limit');axes[1].axhline(.005,color='#777',ls=':')
    axes[0].set(ylabel='Maximum marker position difference');axes[1].set(ylabel='Maximum relative state RMS difference')
    for ax in axes:ax.set_xlabel('Case rank (sorted separately for each comparison)');ax.grid(alpha=.2)
    axes[0].legend();fig.suptitle('Step accuracy and fluid-grid sensitivity are different checks',fontsize=15);save(fig,'05-sensitivity')
    pairs=M['step64']['pairs'];grid=M['grid128']['pairs'];mixed=[p for p in pairs if p['geometry']=='mixed' and p['profile']=='mutual_amplification'];mx={p['initial_state']:p for p in mixed}
    tables=''.join('<tr><td>'+html.escape(p['geometry'])+'</td><td>'+html.escape(p['profile'].replace('_',' '))+'</td><td>'+p['initial_state']+f'</td><td>{p["radius_change_percent"]:+.3f}%</td><td>{q["radius_change_percent"]:+.3f}%</td><td>{p["position_effect_rms_final"]:.6g}</td></tr>' for p,q in zip(pairs,grid))
    geomrows=''.join(f'<tr><td>{n}</td><td>{g["count"]}</td><td>{g["components"]}</td><td>{g["degree_range"]}</td><td>{g["minimum_separation"]:.5g}</td></tr>' for n,g in A['geometry'].items())
    sensrows=''.join(f'<tr><td>{label}</td><td>{s["max_position_difference"]:.6g}</td><td>{s["max_relative_state_rms_difference"]:.6g}</td><td>{s["radius_effect_sign_agreement"]}/45</td><td>{s["max_radius_effect_difference_percentage_points"]:.4g}</td></tr>' for label,s in S.items())
    fig=lambda n,cap:f'<figure><a href="{n}.svg"><img src="{n}.png" alt="{html.escape(cap)}"></a><figcaption>{cap}</figcaption></figure>'
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Menger, Koch and Romeo–Juliet: combined tests</title>
<style>body{{margin:0;background:#fafaf7;color:#172b35;font:17px/1.6 system-ui,sans-serif}}main{{max-width:1160px;margin:auto;padding:40px 25px 90px}}h1{{font-size:42px;line-height:1.12;max-width:900px}}h2{{margin-top:48px}}h3{{margin-top:32px}}.kicker{{font-size:13px;letter-spacing:.1em;text-transform:uppercase;color:#126a76}}.lead{{font-size:22px;max-width:900px}}.box{{border-left:5px solid #d18b22;background:#f1eadb;padding:18px 25px}}.cards{{display:flex;gap:20px;flex-wrap:wrap}}.card{{background:#e8efed;padding:18px 24px;border-radius:10px;min-width:170px}}.card b{{display:block;font-size:30px}}img{{width:100%;height:auto}}figure{{margin:30px 0}}figcaption{{font-size:14px;color:#52626a}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:9px;border-bottom:1px solid #cbd5d7;text-align:left}}th{{background:#e8efed}}.scroll{{overflow-x:auto}}pre{{background:#e8efed;padding:18px;overflow:auto;font-size:15px}}a{{color:#126a76}}details{{margin-top:25px}}code{{font-size:.93em}}footer{{margin-top:50px;font-size:14px;color:#52626a}}</style><main>
<div class="kicker">Local exploratory study · completed through t = 0.10</div><h1>Menger cheese + Koch snowflake + changing relationships</h1>
<p class="lead">The combination is implementable. In this selected interaction model, changing relationship states changes how a marker cloud contracts, expands and rearranges. These results do not establish an improvement to Navier–Stokes or a model of water molecules.</p>
<div class="cards"><div class="card"><b>360</b>case trajectories across four configurations</div><div class="card"><b>192</b>markers in each cloud</div><div class="card"><b>{S['step']['passes']}/90</b>cases pass the step-size screen</div><div class="card"><b>0</b>new fluid simulations</div></div>
<h2>What we combined</h2><p>Five marker arrangements use the same count, centroid and starting RMS radius. Each is tested with three Romeo–Juliet response profiles, three starting state conditions, and motion coupling either off or on. This gives 90 cases and 45 matched comparisons. All cases are repeated with a smaller tracking step, 128-grid saved flow, and coarser velocity snapshots.</p>
{fig('01-geometries','Finite marker arrangements. The Menger cloud samples 192 of 400 depth-2 retained cells; Koch uses all 192 depth-3 vertices. Empty regions are not physical walls. Geometry, embedding and graph structure differ, so this does not isolate fractal dimension.')}
<h2>What happened</h2><p>In the mixed Menger–Koch cloud with mutually amplifying responses, positive starting states changed final cloud radius by <strong>{mx['positive']['radius_change_percent']:+.2f}%</strong> relative to passive markers. Reversing those states changed it by <strong>{mx['negative']['radius_change_percent']:+.2f}%</strong>. Balanced starting states gave <strong>{mx['balanced']['radius_change_percent']:+.2f}%</strong>. These signs describe the measured cloud-size effect, which also includes the spatially varying background flow.</p>
{fig('02-radius-effects','All 45 matched comparisons at both fluid grids. Values compare each active cloud to its own passive counterpart. The three initial states are distinct test conditions, not independent replicates.')}
{fig('03-mixed-motion','Illustrative paths from the saved 64-grid flow. Only every sixth path is drawn; every final marker is shown. Projection crossings do not imply collisions in three dimensions.')}
{fig('04-state-response','Fading, amplification and opposed responses are prescribed by coefficients. An opposing two-person pair has closed linear orbits, but that does not make an entire heterogeneous network a center.')}
<h2>How the equations connect</h2><p>Each marker has position <i>xᵢ</i> and a signed relationship state <i>zᵢ</i>. The initial neighbor graph is symmetric; its weights change with separation. In a two-marker graph, the state equation reduces exactly to the Romeo–Juliet linear system.</p>
<pre>zᵢ′ = aᵢ zᵢ + bᵢ Σⱼ Pᵢⱼ(x) zⱼ
xᵢ′ = u_saved(xᵢ,t) + (μ / d̄₀) Σⱼ wᵢⱼ tanh((zᵢ+zⱼ)/2) rᵢⱼ

rᵢⱼ = shortest periodic displacement from i toward j
wᵢⱼ = exp[−min(|rᵢⱼ|² / (2·0.1²), 50)] on initial graph edges
Pᵢⱼ = wᵢⱼ / Σⱼ wᵢⱼ; d̄₀ = initial mean weighted degree
μ = 0 (passive), or 20 (added attraction/repulsion)</pre>
<p>Positive pair states add attraction; negative pair states add repulsion. The sum of the added pair drifts is zero to rounding accuracy. This is an overdamped drift law: it has no separate particle inertia, collision exclusion or feedback force on the fluid. Zero total drift is not a proof of momentum conservation or incompressibility.</p>
<div class="scroll"><table><tr><th>Profile</th><th>a</th><th>b</th><th>Isolated pair</th></tr><tr><td>Fading</td><td>−20</td><td>10</td><td>Eigenvalues −10, −30: stable node</td></tr><tr><td>Mutual amplification</td><td>−20</td><td>30</td><td>Eigenvalues +10, −50: saddle</td></tr><tr><td>Opposing responses</td><td>0</td><td>+20 and −20</td><td>Opposite pair eigenvalues ±20i: center</td></tr></table></div>
<p>For the full network, half the opposed-response coefficients are +20 and half −20, assigned by a fixed seed. The network geometry matters. The coefficients, coupling strength, connection length and initial states are illustrative choices; none were fitted to water.</p>
<h2>Does this help Navier–Stokes?</h2><div class="box"><strong>Useful as a diagnostic and a separate interaction experiment; no demonstrated improvement to NS.</strong><p>Passive markers can measure how the existing velocity field moves structured neighborhoods. When relationships change marker motion, the added law creates a different particle model. Improving a physical fluid model would require a justified coupling law, compatible conservation and stress relations, convergence, and validation against independent fluid observations. This pilot has no such validation.</p></div>
<p>The effect also occurs in smooth control clouds. A response to the added law therefore cannot by itself establish a special benefit from fractal geometry. The mixed case is a combination of two point arrangements; it is not a new limiting fractal construction.</p>
<h2>Checks, failures and sensitivity</h2><p>All 360 trajectories completed with finite saved values. All 90 cases pass the declared tracking-step screen: largest position difference <strong>{S['step']['max_position_difference']:.3g}</strong>, largest relative state RMS difference <strong>{S['step']['max_relative_state_rms_difference']:.3g}</strong>. The screen limits are 0.001 and 0.005 respectively. Pair integrations agree with matrix-exponential solutions; label permutations, zero-state controls, passive position identity and state-sign symmetry pass. Saved source, input-field and result hashes pass.</p>
{fig('05-sensitivity','Maximum differences over all common saved times. The grid comparison changes the sampled initial fluid as well as spatial resolution; it is sensitivity evidence, not a clean order-of-convergence estimate.')}
<div class="scroll"><table><tr><th>Comparison</th><th>Largest position difference</th><th>Largest relative state difference</th><th>Radius-effect sign agreement</th><th>Largest radius-effect shift (percentage points)</th></tr>{sensrows}</table></div>
<p>The radius-effect sign agrees across grids in 37 of 45 comparisons. All eight reversals occur in balanced-state cases with effects below 0.024% in magnitude. The larger positive/negative-state effects keep their signs. Individual paths nevertheless differ by as much as 0.20 between grids, and by 0.006 when velocity snapshots are coarsened. The mixed cloud's amplification results at 128 grid are −4.72% and +4.94%, compared with −5.11% and +5.24% at 64 grid.</p>
<p><strong>Underlying fluid limitation:</strong> the original raw strain is not smooth across periodic joins, its initial maxima lie outside the central tube, and its strongest gradients are unresolved. Those fields were preserved. Smaller marker-integration errors cannot repair that input limitation. No physical accuracy or asymptotic stability claim follows from passing the tracking screen.</p>
<p>Cloud radius is computed from unwrapped marker coordinates. Separation uses the shortest periodic displacement. Neighbor retention uses the initial directed eight-nearest-neighbor sets, including distance ties; the evolution graph uses their symmetric union. Minimum sampled separation does not prove collision-free motion between saved times. Finite marker clouds do not establish a limiting fractal dimension, and this short window does not estimate asymptotic Lyapunov exponents.</p>
<details><summary>All matched endpoint results</summary><div class="scroll"><table><tr><th>Geometry</th><th>Response</th><th>Initial states</th><th>Radius change, 64</th><th>Radius change, 128</th><th>Position-effect RMS, 64</th></tr>{tables}</table></div></details>
<details><summary>Starting graph checks</summary><table><tr><th>Geometry</th><th>Markers</th><th>Connected components</th><th>Degree range</th><th>Minimum separation</th></tr>{geomrows}</table></details>
<h2>Reproducible files</h2><p><a href="protocol.json">Protocol</a> · <a href="study.py">Simulation script</a> · <a href="analyze.py">Analysis script</a> · <a href="report.py">Report script</a> · <a href="results.json">Run records and hashes</a> · <a href="analysis.json">All measurements and checks</a> · <a href="verification.json">Independent artifact verification</a></p>
<p><a href="base64.npz">Base trajectories</a> · <a href="step64.npz">Smaller-step trajectories</a> · <a href="grid128.npz">128-grid trajectories</a> · <a href="snapshot_coarse.npz">Coarser-snapshot trajectories</a> · <a href="construction.npz">Initial geometry and graph</a></p>
<p>Run analysis and report scripts with <code>work/study-env/Scripts/python.exe</code>. The simulation script requires the existing matched-stretch source and saved field directory, with hashes listed in <code>results.json</code>. Its existing-state guard prevents accidental duplicate runs. Python dependencies: NumPy, SciPy and Matplotlib.</p>
<footer>Exploratory study. All outputs use nondimensional model units. Mathematical model under examination and testing. No added fluid force or new fluid evolution. Publication of this report was subsequently authorized by the user; the original simulation protocol retains its earlier local-only scope as a historical record.</footer></main></html>'''
    (ROOT/'report.html').write_text(page,encoding='utf-8')
    print('Report and five PNG/SVG figure pairs written.')

if __name__=='__main__':main()
