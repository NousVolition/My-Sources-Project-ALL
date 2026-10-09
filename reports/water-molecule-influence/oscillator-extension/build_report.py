"""Build original scientific figures and reports from saved test outputs."""
from pathlib import Path
from html import escape
import base64, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;D=HERE/'data'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.titleweight':'bold','figure.facecolor':'white','axes.labelcolor':'#17313b',
                     'axes.prop_cycle':plt.cycler(color=['#087e83','#df7636','#6258a2','#599749'])})


def main():
    r=json.loads((D/'results.json').read_text())
    e=json.loads((D/'exercises.json').read_text())
    b=json.loads((D/'bifurcations.json').read_text())
    assert r['passed'] and e['passed'] and b['passed']
    z=np.load(D/'trajectories.npz');q=np.load(D/'exercises.npz');h=np.load(D/'bifurcations.npz')
    for archive in [z,q,h]:
        for key in archive.files:assert np.isfinite(archive[key]).all(),key
    def save(fig,name):
        fig.tight_layout(pad=2.)
        fig.savefig(HERE/(name+'.png'),dpi=150,bbox_inches='tight')
        fig.savefig(HERE/(name+'.svg'),bbox_inches='tight')
        plt.close(fig)

    fig,ax=plt.subplots(2,2,figsize=(12,9))
    for radius in [.2,1.,3.]:
        key=f'weak_0.1_{radius:g}'
        t=z[key+'_peak_t'];p=z[key+'_peaks']
        line=ax[0,0].plot(.1*t,p,'o',ms=3,label=f'initial amplitude {radius:g}')[0]
        tau=np.linspace(0,16,501)
        ax[0,0].plot(tau,2/np.sqrt(1+(4/radius**2-1)*np.exp(-tau)),color=line.get_color(),lw=1)
        path=z[key+'_z']
        ax[0,1].plot(path[:,0],path[:,1],lw=.55,alpha=.8)
    ax[0,0].set(title='Van der Pol selects amplitude near 2',xlabel='Slow time μt (μ = 0.1)',ylabel='Positive peak x')
    ax[0,0].legend(fontsize=9);ax[0,0].text(.04,.1,'dots: numerical peaks; lines: averaging',transform=ax[0,0].transAxes,fontsize=9)
    ax[0,1].set(title='Starts inside and outside approach a cycle',xlabel='x',ylabel='v = dx/dt');ax[0,1].set_aspect('equal',adjustable='box')
    for amp in [.5,1.,2.]:
        path=z[f'duffing_{amp:g}_z'];t=z[f'duffing_{amp:g}_t'];mask=t<12
        ax[1,0].plot(path[mask,0],path[mask,1],label=f'amplitude {amp:g}')
    ax[1,0].set(title='Conservative Duffing retains distinct orbits',xlabel='x',ylabel='v');ax[1,0].legend(fontsize=9)
    amps=np.array([row['amplitude'] for row in e['pendulum']]);err=np.array([row['averaging_error'] for row in e['pendulum']])
    ax[1,1].loglog(amps,err,'o-',label='|1 − a²/16 − exact frequency|')
    ax[1,1].loglog(amps,err[0]*(amps/amps[0])**4,'--',label='a⁴ reference')
    ax[1,1].set(title='Small-pendulum frequency error is fourth order',xlabel='Initial angle a (radians)',ylabel='Absolute frequency error');ax[1,1].legend(fontsize=9)
    save(fig,'weak_oscillators')

    fig,ax=plt.subplots(2,2,figsize=(12,9))
    path=z['relax_10_z'];t=z['relax_10_t']
    xx=np.linspace(-2.25,2.25,600)
    ax[0,0].plot(xx,xx**3/3-xx,'k--',lw=1,label='cubic nullcline y = x³/3 − x')
    ax[0,0].plot(path[:,0],path[:,1],lw=1.8,label='μ = 10 cycle')
    for idx in [500,2500,5000,8000,10000]:
        d=path[min(idx+25,len(path)-1)]-path[idx]
        if np.linalg.norm(d)>0:
            unit=d/np.linalg.norm(d)*.13
            ax[0,0].annotate('',xy=path[idx]+unit,xytext=path[idx],arrowprops={'arrowstyle':'->','color':'#df7636'})
    ax[0,0].set(title='Slow branches connected by fast jumps',xlabel='x',ylabel='Transformed coordinate y');ax[0,0].legend(fontsize=8)
    for mu in [10.,40.]:
        ax[0,1].plot(z[f'relax_{mu:g}_t']/mu,z[f'relax_{mu:g}_z'][:,0],label=f'μ = {mu:g}')
    ax[0,1].set(title='One cycle, plotted on the slow time scale',xlabel='Time / μ',ylabel='x');ax[0,1].legend()
    mu=np.array([a['mu'] for a in r['relaxation']])
    slow=np.array([a['slow_interval_duration'] for a in r['relaxation']]);fast=np.array([a['fast_interval_duration'] for a in r['relaxation']])
    ax[1,0].loglog(mu,slow,'o-',label='slow: x 1.8 → 1.2')
    ax[1,0].loglog(mu,fast,'o-',label='fast: x 0.5 → −0.5')
    ax[1,0].loglog(mu,r['asymptotic_constants']['slow_over_mu']*mu,'k--',lw=.8,label='asymptotic ∝ μ and 1/μ')
    ax[1,0].loglog(mu,r['asymptotic_constants']['fast_times_mu']/mu,'k--',lw=.8)
    ax[1,0].set(title='Two explicitly defined traversal times',xlabel='μ',ylabel='Duration');ax[1,0].legend(fontsize=8)
    ax[1,1].plot(mu,[a['period']/a['mu'] for a in r['relaxation']],'o-',label='numerical T/μ')
    ax[1,1].axhline(r['asymptotic_constants']['period_over_mu'],color='black',ls='--',label='3 − 2 ln 2')
    ax[1,1].set(title='Period approaches the large-μ prediction',xlabel='μ',ylabel='Period / μ');ax[1,1].legend(fontsize=9)
    save(fig,'relaxation')

    fig,ax=plt.subplots(2,2,figsize=(12,9))
    t=q['cubic_damping_2_t'];path=q['cubic_damping_2_z'];approx=q['cubic_damping_2_approximation']
    ax[0,0].plot(t,path[:,0],label='full cubic velocity damping',lw=1.3)
    ax[0,0].plot(t,approx,'--',label='averaged amplitude × cos t',lw=1)
    ax[0,0].set(title='Supplied challenge: ε = 2, a = 1',xlabel='Time',ylabel='x');ax[0,0].legend(fontsize=8)
    ax[0,1].plot(t,path[:,0]-approx)
    ax[0,1].set(title='Approximation error is measured, not assumed',xlabel='Time',ylabel='Full x − averaged x')
    for seed in [0.,1e-4]:
        t=q[f'swing_{seed:g}_t'];path=q[f'swing_{seed:g}_z']
        ax[1,0].plot(t,np.linalg.norm(path,axis=1),label=f'initial x = {seed:g}, v = 0')
    ax[1,0].set(title='Pumping amplifies a seed; exact rest stays at rest',xlabel='Time',ylabel='√(x² + v²)');ax[1,0].legend(fontsize=8)
    for eps in [.1,.05,.025]:
        rows=[row for row in e['swing']['floquet'] if row['epsilon']==eps]
        ax[1,1].plot([row['gamma'] for row in rows],[row['growth_rate_per_time']/eps for row in rows],'o',label=f'ε = {eps:g}')
    gamma=np.linspace(-.8,.8,501)
    ax[1,1].plot(gamma,np.sqrt(np.maximum(0,1-4*gamma*gamma))/4,'k--',label='leading averaging prediction')
    ax[1,1].set(title='Growth band from the one-period linear map',xlabel='Detuning γ',ylabel='Exponential growth rate / ε');ax[1,1].legend(fontsize=8)
    save(fig,'exercise_controls')

    fig,ax=plt.subplots(1,2,figsize=(12,4.5))
    mu=np.linspace(-.3,.3,601);pos=mu>=0
    ax[0].plot(mu[~pos],np.zeros(np.sum(~pos)),label='stable equilibrium')
    ax[0].plot(mu[pos],np.zeros(np.sum(pos)),'--',color='#df7636',label='unstable equilibrium')
    ax[0].plot(mu[pos],np.sqrt(mu[pos]),color='#087e83');ax[0].plot(mu[pos],-np.sqrt(mu[pos]),color='#087e83')
    ax[0].set(title='Pitchfork: an eigenvalue is zero at μ = 0',xlabel='μ',ylabel='Equilibrium x');ax[0].legend(fontsize=9)
    for val in [-.2,0.,.2]:
        key=f'hopf_{val:g}_0.1';path=h[key+'_z']
        ax[1].plot(h[key+'_t'],np.linalg.norm(path,axis=1),label=f'μ = {val:g}')
    ax[1].axhline(np.sqrt(.2),color='black',ls='--',lw=.8,label='cycle radius √0.2')
    ax[1].set(title='Hopf: eigenvalues cross at ±i, determinant = 1',xlabel='Time',ylabel='Radius from initial r = 0.1');ax[1].legend(fontsize=9)
    save(fig,'bifurcations')

    rr=r['relaxation'][-1];damping=e['cubic_velocity_damping']['rows'][-1]
    swing=e['swing']['nonlinear_seed_test'][-1]
    weak01=[row for row in r['weak'] if row['mu']==.1]
    findings=[
       ['Weak van der Pol, μ=0.1',f'Three initial amplitudes (0.2, 1, 3) give last positive peaks {min(a["last_peak"] for a in weak01):.6f}–{max(a["last_peak"] for a in weak01):.6f} by μt≈16.'],
       ['Relaxation at μ=40',f'Slow interval {rr["slow_interval_duration"]:.6f}; central jump {rr["fast_interval_duration"]:.6f}; ratio {rr["slow_fast_ratio"]:.2f}.'],
       ['Large-μ period',f'T={rr["period"]:.6f}, compared with μ(3−2 ln 2)={rr["period_asymptote"]:.6f}; {100*rr["period_relative_asymptotic_error"]:.3f}% asymptotic error.'],
       ['Conservative Duffing',f'Three distinct amplitudes persist. Largest relative energy error over t≤200: {max(a["maximum_relative_energy_error"] for a in r["duffing"]):.2e}.'],
       ['Cubic velocity damping, ε=2',f'Over 0≤t≤50, averaged-waveform maximum error {damping["maximum_position_error"]:.6f}, RMS error {damping["RMS_position_error"]:.6f}.'],
       ['Pumped swing',f'Exact rest stays zero. A 0.0001 seed reaches radius {swing["final_radius"]:.6f} at t=160, ε=0.1, γ=0.'],
       ['Numerical controls',f'RK4 observed orders {r["rk4"]["orders"][0]:.4f}, {r["rk4"]["orders"][1]:.4f}; stiff-solver tolerance refinement and a second solver agree.']]
    relaxation_table=[[f'{a["mu"]:g}',f'{a["period"]:.6f}',f'{100*a["period_relative_asymptotic_error"]:.3f}%',
                      f'{a["slow_interval_duration"]:.6f}',f'{a["fast_interval_duration"]:.6f}',f'{a["slow_fast_ratio"]:.2f}'] for a in r['relaxation']]
    sections=[
      ('Question and scope','Can the same van der Pol equation show weak amplitude selection and strong fast–slow relaxation, and do the controls separate a stable limit cycle, conserved closed orbits, damping, parametric growth, and a Hopf bifurcation? These are dimensionless mathematical benchmarks. They add no molecular-dynamics measurements and establish no water-molecule agency or hierarchy.'),
      ('Equations and coordinate conventions','The supplied van der Pol equation is x″+μ(x²−1)x′+x=0. For weak μ we integrate (x,v), with v=x′. For large μ, set F(x)=x³/3−x and y=F(x)+v/μ; then x′=μ[y−F(x)] and y′=−x/μ. Thus the plotted y is not velocity. The relaxation tests begin at (x,y)=(2,0), corresponding to v=−2μ/3; the weak tests begin at v=0.'),
      ('Weak van der Pol: a testable amplitude prediction','Averaging gives dr/dt=(μr/2)(1−r²/4), so r(t)=2/[1+(4/r₀²−1)e^(−μt)]^(1/2). We compare this envelope with event-located positive maxima, not with the rapidly oscillating signed x. Use μ=0.05, 0.1, 0.2; r₀=0.2, 1, 3; run until μt=16. At μ=0.05 the largest peak-envelope discrepancy is '+f'{max(a["max_peak_envelope_error"] for a in r["weak"] if a["mu"]==.05):.6f}'+'. Finite-μ cycles are not exactly circles, and their peaks are not exactly 2.'),
      ('Green’s-theorem check of the approximate radius','In (x,v), the van der Pol field has divergence μ(1−x²). A periodic orbit has zero net outward flow through its boundary. Approximating its enclosed region by a disk of radius r gives integral μ π r²(1−r²/4), whose nonzero root is r=2. The numerical area integrals for r=1,2,3 agree with this expression. The disk assumption is approximate; zero flux alone does not prove existence or attraction of a periodic orbit.'),
      ('What is fast, and what is slow?','The outer cubic branches |x|>1 attract fast horizontal motion because the fast x derivative is μ(1−x²)<0. On the slow portion, x′≈−x/[μ(x²−1)], so traversal takes O(μ). Away from the folds, y−F(x) is O(μ⁻²), not exactly zero. During the central part of a jump, x′ is O(μ), giving an O(μ⁻¹) traversal time. Near a fold these simple approximations need a boundary-layer analysis.'),
      ('Precisely defined time-scale measurements','Within the last converged cycle, measure the descending outer-branch segment x=1.8→1.2 and the descending central jump x=0.5→−0.5. These are portions of the motion, not the complete slow leg or the complete jump including departure from the fold. Event times come from bracketed roots of the dense solution. The asymptotic constants are slow time/μ=0.9−ln(1.5)='+f'{r["asymptotic_constants"]["slow_over_mu"]:.6f}'+', and μ×fast time='+f'{r["asymptotic_constants"]["fast_times_mu"]:.6f}'+'. At μ=40 the measured values are '+f'{rr["slow_over_mu"]:.6f} and {rr["fast_times_mu"]:.6f}'+'. The median slow-branch nullcline offset falls from '+f'{r["relaxation"][0]["median_slow_nullcline_residual"]:.6f} at μ=5 to {rr["median_slow_nullcline_residual"]:.8f} at μ=40.'),
      ('Relaxation period and solver controls','Integrating the singular-limit slow branches gives T≈μ(3−2 ln 2). Its error decreases from '+f'{100*r["relaxation"][0]["period_relative_asymptotic_error"]:.2f}% at μ=5 to {100*rr["period_relative_asymptotic_error"]:.2f}% at μ=40.'+' Each case is solved with Radau at relative tolerances 10⁻⁹ and 10⁻¹¹; periods and consecutive late cycles must agree within 10⁻⁶ relative. An independent DOP853 solution at μ=10 agrees along the sampled final cycle within '+f'{r["relaxation"][1]["DOP853_cycle_max_distance"]:.2e}'+'. This separates asymptotic-model error from observed solver disagreement. It does not certify every value to all printed digits.'),
      ('Duffing control: periodic does not mean attracting','For the supplied x″+x+εx³=0, E=v²/2+x²/2+εx⁴/4 is conserved. With ε=0.1 and initial amplitudes 0.5,1,2, the trajectories retain their distinct amplitudes through t=200. Positive-peak periods agree with an independent energy quadrature to '+f'{max(a["period_relative_error"] for a in r["duffing"]):.2e}'+ ' relative. First-order averaging predicts frequency 1+3εa²/8; its error increases with εa². These nested closed orbits are not attracting limit cycles.'),
      ('Pendulum frequency control','For x″+sin x=0, the leading small-angle frequency is 1−a²/16. Test a=0.1,0.2,0.4 against the exact period 4K(sin²(a/2)), where K uses the elliptic parameter convention. The frequency errors are '+', '.join(f'{a["averaging_error"]:.3e}' for a in e['pendulum'])+'. Their measured amplitude-scaling orders are '+', '.join(f'{x:.3f}' for x in e['pendulum_averaging_error_orders'])+', consistent with the omitted fourth-order term.'),
      ('The ε=2 exercise is cubic velocity damping','The printed equation is x″+ε(x′)³+x=0, not the conservative Duffing equation. With x(0)=a, v(0)=0, first-order averaging gives r(t)=a/[1+3εa²t/4]^(1/2), and x(t)≈r(t)cos t. We test a=1, ε=0.1,0.2,2 through t=50. At ε=2 the initial energy 0.5 falls to '+f'{damping["final_energy"]:.6f}'+'. The independently accumulated loss ∫εv⁴dt closes the energy balance to '+f'{damping["maximum_energy_balance_residual"]:.2e}'+'. The waveform comparison is plotted and its actual errors are reported; ε=2 is outside a controlled small-ε expansion. The leading approximation also has an O(ε) initial-velocity mismatch from its decaying envelope.'),
      ('Swing: a push and an instability are different','For x″+[1+εγ+ε cos(2t)]sin x=0, x=v=0 is an exact solution. Linearizing at rest and averaging with slow time T=εt gives dr/dT=r sin(2φ)/4 and dφ/dT=[γ+cos(2φ)/2]/2. The leading instability band is |γ|<1/2, with physical-time growth ε sqrt(1−4γ²)/4. We integrate the 2×2 fundamental matrix over the forcing period π and measure the eigenvalue growth. At ε=0.1, γ=0 the numerical exponent is '+f'{next(a["growth_rate_per_time"] for a in e["swing"]["floquet"] if a["epsilon"]==.1 and a["gamma"]==0):.8f}'+', compared with 0.025 from averaging. Cases γ=±0.75 remain outside the sampled growth band. Finite-ε band boundaries can shift; these sampled points do not locate the exact boundaries. In the full nonlinear model, a tiny seed grows while exact rest remains zero. Real disturbances supply seeds, but no noise was added here.'),
      ('Zero eigenvalue versus Hopf bifurcation','The latest cropped phase portrait does not specify a complete equation, so we use two explicitly chosen normal forms. Pitchfork: x′=μx−x³, y′=−y, with eigenvalues (μ,−1) at the origin. Hopf: x′=μx−y−xr², y′=x+μy−yr², r²=x²+y², with eigenvalues μ±i. At the Hopf threshold μ=0, determinant=1: neither eigenvalue is zero. In polar coordinates r′=μr−r³, θ′=1. For μ>0 the stable cycle has radius sqrt(μ) and period 2π. Numerical trajectories for μ=−0.2,0,0.2 and initial radii 0.1,1 agree with exact radial solutions within '+f'{max(a["exact_max_distance"] for a in b["hopf"]["rows"]):.2e}'+'. At μ=0 the cubic term still makes the origin asymptotically stable. This is a supercritical example, not a claim that every Hopf bifurcation creates a stable cycle.'),
      ('Limitations and next validation','These finite integrations test stated equations, chosen starting states and explicit numerical controls. They do not prove global behavior for arbitrary models. Next tests: estimate a return-map multiplier for van der Pol attraction; resolve fold departure with larger μ; locate finite-ε swing instability boundaries; add controlled noise and quantify phase diffusion; and supply the full equation for any cropped diagram to test that particular system. Linking any oscillator model to water would require a separately specified physical mechanism and molecular validation.')]
    design=[['Amplitude selection','3 μ values × 3 initial amplitudes; positive-peak events through μt=16.','Compare analytic averaged envelope; last peaks within 0.01 of 2; spread below 0.002.'],
            ['Fast–slow separation','μ=5,10,20,40; final full cycle after t up to 12μ+40.','Refine Radau tolerance; compare DOP853 at μ=10; use fixed traversal thresholds.'],
            ['Conservation and damping','Duffing ε=0.1; cubic velocity damping ε=0.1,0.2,2.','Energy quadrature or accumulated loss; no assumed agreement outside small ε.'],
            ['Small-angle prediction','Pendulum a=0.1,0.2,0.4.','Exact elliptic period; omitted frequency term scales as a⁴.'],
            ['Parametric growth','Swing ε=0.1,0.05,0.025; γ=−0.75,−0.25,0,0.25,0.75.','One-period linear map; determinant=1; nonlinear zero/seed controls.'],
            ['Bifurcation mechanism','Explicit pitchfork and supercritical Hopf normal forms.','Eigenvalues, determinant, exact radius evolution; do not infer missing equations.']]
    def mdtable(head,rows):
        clean=lambda v:str(v).replace('|','\\|')
        return '| '+' | '.join(head)+' |\n| '+' | '.join(['---']*len(head))+' |\n'+'\n'.join('| '+' | '.join(clean(v) for v in row)+' |' for row in rows)
    def table(head,rows):
        return '<div class="scroll"><table><tr>'+''.join('<th>'+escape(x)+'</th>' for x in head)+'</tr>'+''.join('<tr>'+''.join('<td>'+escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</table></div>'
    commands='''python -m pip install -r oscillator-extension/requirements.txt
python oscillator-extension/test_oscillators.py
python oscillator-extension/test_pendulum_swing.py
python oscillator-extension/test_bifurcations.py
python oscillator-extension/build_report.py'''
    title='Oscillators, damping, and bifurcations: controlled numerical tests'
    readme='# '+title+'\n\n'+mdtable(['Test','Recorded result'],findings)+'\n\n## Test design\n\n'+mdtable(['Stage','Measurement','Control'],design)
    placements={2:'weak_oscillators',6:'relaxation',10:'exercise_controls',11:'bifurcations'}
    for idx,(heading,body) in enumerate(sections):
        readme+='\n\n## '+heading+'\n\n'+body
        if idx in placements:readme+='\n\n!['+heading+']('+placements[idx]+'.png)'
    readme+='\n\n## Relaxation measurements\n\n'+mdtable(['μ','Period','Asymptotic error','Slow time','Fast time','Ratio'],relaxation_table)
    readme+='\n\n## Reproduce\n\nUse Python 3.12, from the parent experiment folder:\n\n```sh\n'+commands+'\n```\n\nNo random seed is needed. The scripts save numeric JSON and compressed trajectory arrays, stop on failed numerical acceptance checks, and rebuild these original figures and reports. Small differences across numerical-library versions are expected.\n\n[Full report](report.html) · [Protocol](data/protocol.json) · [Oscillator results](data/results.json) · [Exercise results](data/exercises.json) · [Bifurcation results](data/bifurcations.json) · [Checksums](manifest_sha256.json)\n\nEquations are transcribed from the supplied excerpts, except the clearly labeled chosen normal forms. The scanned book pages are not redistributed. [Earlier stability tests](../stability-extension/README.md) · [Original water study](../README.md)\n'
    (HERE/'README.md').write_text(readme,encoding='utf-8')
    html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Controlled oscillator tests</title><style>body{font:17px/1.65 system-ui,sans-serif;max-width:1120px;margin:40px auto;padding:0 24px;background:#fafcfd;color:#17313b}h1{font-size:40px;line-height:1.2}h2{margin-top:38px;font-size:25px}.box{border-left:4px solid #087e83;background:#e8f3f3;padding:18px}.scroll{overflow-x:auto}table{width:100%;border-collapse:collapse;font-size:14px}th,td{padding:10px;text-align:left;border-bottom:1px solid #d7e3e6}th{background:#e8f0f2}img{max-width:100%;height:auto;margin:20px 0}pre{background:#edf3f5;padding:20px;overflow:auto}a{color:#087e83}</style></head><body><p>REPRODUCIBLE MATHEMATICAL BENCHMARKS · 9 OCTOBER 2026</p><h1>Slow motion, fast jumps,<br>and the birth of oscillations</h1><p class="box">Different mechanisms can produce similar-looking motion. Test whether amplitude is selected, energy is conserved or lost, a small seed grows, and the prediction survives numerical refinement.</p>'''
    html+=table(['Test','Recorded result'],findings)+'<h2>Test design</h2>'+table(['Stage','Measurement','Control'],design)
    for idx,(heading,body) in enumerate(sections):
        html+='<h2>'+escape(heading)+'</h2><p>'+escape(body)+'</p>'
        if idx in placements:
            name=placements[idx];encoded=base64.b64encode((HERE/(name+'.png')).read_bytes()).decode()
            html+='<img alt="'+escape(heading)+'" src="data:image/png;base64,'+encoded+'">'
    html+='<h2>Relaxation measurements</h2>'+table(['μ','Period','Asymptotic error','Slow time','Fast time','Ratio'],relaxation_table)
    html+='<h2>Reproduce</h2><p>Python 3.12; run from the parent experiment folder. Parameters, raw trajectories, numeric results, and checksums accompany this report.</p><pre>'+escape(commands)+'</pre><p><a href="README.md">Methods and file index</a> · <a href="../stability-extension/report.html">Earlier stability report</a></p><p>Original figures; supplied scanned pages are not redistributed.</p></body></html>'
    (HERE/'report.html').write_text(html,encoding='utf-8')
    (HERE/'TEST_DESIGN.md').write_text('# Test design\n\n'+mdtable(['Stage','Measurement','Control'],design)+'\n\nSee data/protocol.json and the saved exercise/bifurcation results for exact parameters, formulas and acceptance criteria. These are mathematical benchmark tests, not water or biological experiments.\n',encoding='utf-8')
    (D/'audit.json').write_text(json.dumps({'finite_saved_arrays':True,'figure_count':4,'all_three_test_suites_passed':True},indent=2),encoding='utf-8')
    files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest_sha256.json']
    manifest={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    (HERE/'manifest_sha256.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({'files':len(files)+1,'figures':4,'all_checks_passed':True}))


if __name__=='__main__':main()
