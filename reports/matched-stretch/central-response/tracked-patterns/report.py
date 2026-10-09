"""Build reproducible figures and report from saved pair measurements."""
from pathlib import Path
import os,json
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def read(n):
    if n=='measurements.json' and not (ROOT/n).exists():
        manifest=json.loads((ROOT/'data-manifest.json').read_text(encoding='utf-8'))
        data=json.loads((ROOT/manifest['common']).read_text(encoding='utf-8'))
        data['runs']={rid:json.loads((ROOT/path).read_text(encoding='utf-8')) for rid,path in manifest['runs'].items()}
        return data
    return json.loads((ROOT/n).read_text(encoding='utf-8'))
def main():
    data=read('measurements.json');local=read('linearization.json');tracks=read('rk4-tracks.json')
    kinds=('central','outer','across');colors=('#146c94','#bf6b28','#6554a4');comparison=[]
    for case in ('aligned','compressive'):
        for label,a,b in [('grid',f'{case}-fourier-n112-base',f'{case}-fourier-n160-base'),('time step',f'{case}-fourier-n112-base',f'{case}-fourier-n112-half'),('method',f'{case}-fd4-n160-base',f'{case}-fourier-n160-base')]:
            aa=data['runs'][a]['frames'];bb=data['runs'][b]['frames']
            for end in (.1,.4):
                count=2 if end==.1 else 5
                pa=np.array([f['positions'] for f in aa[:count]]);pb=np.array([f['positions'] for f in bb[:count]])
                pos=np.linalg.norm((pa-pb+3)%6-3,axis=2)
                va=np.array([f['velocity'] for f in aa[:count]]);vb=np.array([f['velocity'] for f in bb[:count]])
                comparison.append(dict(case=case,comparison=label,through=end,max_position_difference=float(pos.max()),velocity_relative_l2=float(np.linalg.norm(va-vb)/np.linalg.norm(vb))))
    fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,case in enumerate(('aligned','compressive')):
        frames=data['runs'][f'{case}-fourier-n160-base']['frames'];t=[f['t'] for f in frames]
        for ax,key,title in zip(axs[row],('distance','relative_speed','shared_speed'),('Pair distance','Difference in velocity','Shared speed')):
            for kind,color in zip(kinds,colors):ax.plot(t,[f['summary'][kind][key] for f in frames],'o-',label=kind,color=color,ms=4)
            ax.set(title=f'{case.capitalize()}: {title}',xlabel='Model time',ylabel='Median across fixed pairs');ax.grid(alpha=.2)
            if key=='distance':ax.set_yscale('log')
            ax.legend(fontsize=8)
    fig.suptitle('The same labeled pairs can approach, separate, and change speed | 160 cubed')
    fig.savefig(ROOT/'pair-patterns.png',dpi=155);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4.8),layout='constrained')
    for ax,case in zip(axs,('aligned','compressive')):
        rows=local['runs'][f'{case}-fourier-n160-base']
        for factor,color in zip(('1.0','0.5','0.25'),colors):ax.plot([r['t'] for r in rows],[100*r['probe_comparisons'][factor]['central']['relative_rms_error'] for r in rows],'o-',label=f'{factor} x pair separation',color=color)
        ax.set(title=case.capitalize(),xlabel='Model time',ylabel='Central-pair relative RMS error (%)');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Local matrix prediction vs sampled relative velocity | interpolation included')
    fig.savefig(ROOT/'local-matrix.png',dpi=155);plt.close(fig)
    positions=np.array([f['positions'] for f in tracks['runs']['eighth_step']['frames']]);seed=np.array(tracks['initial_positions'])
    fig,axs=plt.subplots(1,2,figsize=(10,5),layout='constrained')
    for ax,axis,label in ((axs[0],1,'y'),(axs[1],2,'z')):
        for i in range(27):
            group=int(np.argmin(abs(np.array([-.3,0,.3])-seed[i,2])));ax.plot(positions[:,i,0],positions[:,i,axis],color=colors[group],lw=.9,alpha=.8)
            ax.scatter(positions[0,i,0],positions[0,i,axis],s=10,facecolors='none',edgecolors=colors[group]);ax.scatter(positions[-1,i,0],positions[-1,i,axis],s=14,marker='x',color=colors[group])
        ax.set(xlabel='Unwrapped x (domain units)',ylabel=f'Unwrapped {label} (domain units)',title=f'x-{label} position projection');ax.grid(alpha=.2);ax.set_aspect('equal',adjustable='datalim')
    from matplotlib.lines import Line2D
    axs[0].legend(handles=[Line2D([0],[0],color=c,label=f'Initial z = {z}') for c,z in zip(colors,(-.3,0,.3))],fontsize=8,loc='lower right')
    fig.suptitle('27 RK4 tracks in the original saved flow | open circle: start; x: end')
    fig.savefig(ROOT/'rk4-trajectories.png',dpi=155);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    for ax,case in zip(axs,('aligned','compressive')):
        frames=data['runs'][f'{case}-fourier-n160-base']['frames']
        for pair,color in zip(((20,21),(80,81),(140,141)),colors):
            e=next(i for i,e in enumerate(data['edges']) if (e['i'],e['j'])==pair)
            x=[f['metrics']['distance'][e] for f in frames];y=[f['metrics']['radial_speed'][e] for f in frames];ax.plot(x,y,'o-',color=color,label=f'Labels {pair[0]} / {pair[1]}')
            for k in range(4):ax.annotate('',xy=(x[k+1],y[k+1]),xytext=(x[k],y[k]),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
        ax.axhline(0,color='.5',lw=.7);ax.set(title=case.capitalize(),xlabel='Pair distance (domain units)',ylabel='Separation speed (units / model time)');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Observable-plane paths | negative: approaching; positive: separating')
    fig.savefig(ROOT/'pair-phase-plane.png',dpi=155);plt.close(fig)
    rows=[]
    for case in ('aligned','compressive'):
        f=data['runs'][f'{case}-fourier-n160-base']['frames'];rr=local['runs'][f'{case}-fourier-n160-base']
        rows.append({'case':case,'initial_central_pair_distance':f[0]['summary']['central']['distance'],'final_central_pair_distance':f[-1]['summary']['central']['distance'],'final_outer_pair_distance':f[-1]['summary']['outer']['distance'],'final_across_pair_distance':f[-1]['summary']['across']['distance'],'final_field_change':rr[-1]['velocity_field_relative_change_from_start'],'minimum_sampled_label_separation':min(r['minimum_label_separation'] for r in rr)})
    summary={'labels_per_run':192,'pairs_per_run':len(data['edges']),'runs':len(data['runs']),'fields':len(data['verified_fields']),'new_fluid_runs':0,'endpoints':rows,'refinement':comparison,'rk4':tracks['verification']}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    text=['# Relationships tracked in the fluid','','**Neighbors can move in the same direction while their spacing and speeds change. The central labels are already directionally aligned at the start; this is not newly acquired alignment. Later detailed paths remain sensitive to numerical resolution.**','','[Explore labeled pairs](index.html) · [Back to central-tube results](../README.md)','','## What is tracked','','The separate periodic surroundings runs already saved 192 moving labels: 64 along the central axis and 32 along each of four surrounding axes. This analysis follows 377 fixed pairs: 64 central-axis neighbors, 128 outer-axis neighbors, and 185 links between axes. Along each initial axis, consecutive labels are linked with periodic closure. Each label also links to its closest label on another axis at time zero; ties use the lowest label index. Duplicate links are removed. The links remain fixed as the fluid evolves. Across-axis links are relationships between sampled points, not claims that they are molecular neighbors.','','The labels move according to the existing fluid velocity. No pair force, emotional response law or change to Navier–Stokes is introduced. The textbook word “personality” can describe a measured response pattern; the measurements here use distance, direction and speed. A label represents a fluid parcel, not a resolved molecule.','','## The relationship equations','','For particle labels `i,j`, `dX_i/dt = u(X_i,t)` and `d(X_j-X_i)/dt = u(X_j,t)-u(X_i,t)`. Thus relative motion is already determined by the fluid. The distance uses the shortest periodic displacement `r`, and `delta_u = u_j-u_i`. Away from changes of shortest periodic image:','','- Separation speed: `r dot delta_u / |r|`. Negative means approaching; positive means separating.','- Direction agreement: `u_i dot u_j / (|u_i||u_j|)`. +1 is parallel, -1 is opposite. It is undefined when either speed is near zero.','- Difference in velocity: `|u_j-u_i|`.','- Shared speed: `|(u_i+u_j)/2|`.','- Turning rate of the pair separation: `|delta_u - r_hat*(r_hat dot delta_u)| / |r|`.','','Direction agreement and shared speed depend on the chosen reference frame. Pair distance and relative velocity are invariant under adding uniform translation. Smooth rigid rotation can produce opposite particle velocities while maintaining separation; the checks include this case. All quantities are in the existing model units.','','## Measured patterns','','![Pair measurements](pair-patterns.png)','','| Initial surroundings | Initial central spacing | Final central spacing | Final outer spacing | Final across-axis spacing |','| --- | ---: | ---: | ---: | ---: |']
    for r in rows:text.append(f"| {r['case']} | {r['initial_central_pair_distance']:.5f} | {r['final_central_pair_distance']:.5f} | {r['final_outer_pair_distance']:.5f} | {r['final_across_pair_distance']:.5f} |")
    text+=['','Values are medians over the same fixed pairs, not widths of a vortex core or density measurements. Incompressibility permits contraction along one direction with expansion in others. The medians of the central direction cosines are +1 at all saved times; four pairs touching stationary labels have undefined direction comparison. That symmetry was present at initialization. Outer neighbors also show high direction agreement, but both relative and shared speeds change. Medians can hide differently behaving pairs, which the interactive view exposes.','','The full velocity field changes from its initial state by '+', '.join(f"{100*r['final_field_change']:.1f}% ({r['case']})" for r in rows)+' in relative spatial L2 norm at time 0.40. These are evolving flows, not a demonstrated steady pattern.','','## Does the local linear model work?','','![Local linearization check](local-matrix.png)','','For nearby points, `delta_u` is approximated by `G*r`, where `G` is the velocity-gradient matrix at the pair midpoint. Here G is the exact derivative inside a cell of the trilinear interpolation used to read saved velocity; it changes at cell faces. We compare that prediction with the interpolated velocity difference, then repeat using half and quarter separations about the same midpoint. These shorter probes are diagnostic positions, not newly evolved particles.','','The initial central-pair relative RMS discrepancy is about 0.8%. At time 0.40, the aligned case has 27.1% discrepancy at the full separation and 3.2% at quarter separation; the compressive case changes from 9.5% to 10.1%. Smaller neighborhoods do not give a consistent monotone improvement in this test. The discrepancy includes the piecewise interpolation and variation of the velocity field. It is not a measurement of the Navier–Stokes nonlinear term alone.','','No fixed point or autonomous two-variable model was identified. These matrices therefore do not establish stable/unstable manifolds or asymptotic stability. In a smooth autonomous system, a hyperbolic saddle is the setting where the textbook local saddle conclusion applies. The moving, time-dependent fluid comparison requires its own trajectory checks.','','## Trajectories and the phase-plane connection','','![Observable-plane paths](pair-phase-plane.png)','','These are paths through measured distance and separation speed. Arrows indicate time order. They are projections of a larger evolving system, not a closed two-dimensional phase portrait; crossings do not contradict uniqueness.','','![Runge–Kutta tracks](rk4-trajectories.png)','','Separately, 27 new passive tracks were integrated with RK4 through the **original matched-stretch 64³ saved velocity** from time 0 to 0.40. Initial positions are all combinations of `x,y = -0.15,0,0.15` and `z = -0.3,0,0.3`. They use periodic trilinear spatial interpolation and linear interpolation between velocity snapshots spaced by 0.01. The original fluid solution is not recomputed. These start and field differ from the 192-label surroundings runs above.','','The plotted x-y and x-z views are physical-position projections. A crossing can occur at a different depth or time. They are not autonomous planar phase portraits.','','| RK4 step | Half step | Maximum position difference |','| --- | --- | ---: |']
    for r in tracks['verification']['step_ladder']:text.append(f"| {r['dt']} | {r['half_dt']} | {r['max_position_difference']:.7f} |")
    v=tracks['verification'];text += ['',f"At the finest step comparison, the maximum difference is {v['step_halving_max_position_difference']:.7f} domain units (cell size {v['spatial_cell_size']}). Coarsening the saved-velocity interval from 0.01 to 0.02 changes tracks by up to {v['saved_time_coarsening_max_position_difference']:.5f}; this is a sensitivity test, not an error bound for the original snapshots. Through time 0.10 those two differences are {v['through_0_1_step_difference']:.3g} and {v['through_0_1_saved_time_difference']:.5f}. Later plotted tracks must be read as diagnostic reconstructions.",'','## Separation and refinement','','No two of the 192 stored labels occupy exactly the same position at the five checked times. The smallest sampled separations are '+', '.join(f"{r['minimum_sampled_label_separation']:.6f} ({r['case']})" for r in rows)+'. These separations can be below a grid cell: a passive interpolation label can move continuously inside a cell, but this does not mean the underlying velocity variation at that scale is resolved. The snapshots do not prove no collision occurred between outputs or establish uniqueness of a limiting continuum solution.','','| Case | Comparison | Through time | Maximum corresponding-label position difference | Relative sampled-velocity L2 difference |','| --- | --- | ---: | ---: | ---: |']
    for r in comparison:text.append(f"| {r['case']} | {r['comparison']} | {r['through']:.2f} | {r['max_position_difference']:.6f} | {100*r['velocity_relative_l2']:.3f}% |")
    text += ['','Grid compares 112³ with 160³; time step compares 112³ base with half step; method compares FD4 and Fourier at 160³. The two methods share pressure projection/filter infrastructure. All comparisons use the same labels at the same saved times. The existing late peak-resolution failures remain in force. Original matched-stretch also retains its raw periodic-join limitation. Neither a finite set of particle labels nor small time-step error establishes a resolved fluid peak.','','## Verification and reproduction','','40 immutable surroundings fields were hashed, source versions checked, and all 192 velocity samples in each field compared with independent SciPy interpolation. Analytic controls check periodic distance, uniform translation, rigid rotation and linear extension. The separate RK4 control follows a known circular orbit and reduces its error by about 16 when halving the step. Original RK4 tracks read and hash 41 saved velocity fields. No running numerical source was edited.','','The data files contain numerical measurements and new analysis only. Private uploaded reports and screenshots are excluded. Browser preview was unavailable under the existing local-file URL restriction; data and script checks are recorded separately from visual browser verification.','','[Summary](summary.json) · [Local matrix measurements](linearization.json) · [RK4 measurements](rk4-tracks.json) · [Pair-data manifest](data-manifest.json) · [Analysis code](analyze.py) · [Report builder](report.py)','']
    text+=['## Known reference: the double well','','[View the separately tested textbook equation and its phase portrait](double-well-reference/README.md). RK4 checks distinguish closed motion, the saddle boundary and travel around both wells. This supplies a known-answer comparison; no double-well potential is assigned to the fluid.','']
    (ROOT/'README.md').write_text('\n'.join(text),encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
