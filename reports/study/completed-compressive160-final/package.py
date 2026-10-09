"""Package the completed compressive controls from saved data; no simulation."""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, shutil, os, ast, re
import numpy as np
HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
PUB=HERE/'publication'
DEST=PUB/'reports/study'
BATCH=DEST/'completed-compressive160-final'
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
IDS=('compressive-fourier-n160-base','compressive-fourier-n160-half','compressive-fd4-n160-base','compressive-fd4-n160-half')
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def blob(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,s):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s,encoding='utf-8')
def dump(p,d):write(p,json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
    base={v['path']:v for v in read(HERE/'base-tree.json')['tree'] if v['type']=='blob'}
    def base_text(name,remote):
        raw=(HERE/name).read_bytes()
        lf=raw.replace(b'\r\n',b'\n')
        candidates=[s for n in range(3) for s in (lf[:len(lf)-n],lf[:len(lf)-n].replace(b'\n',b'\r\n'))]
        raw=next(b for b in candidates if blob(b)==base[remote]['sha'])
        return raw.decode('utf-8').replace('\r\n','\n')
    old_study=base_text('base-study-README.md','reports/study/README.md')
    old_index=base_text('base-reports-README.md','reports/README.md')
    old_batch=base_text('base-previous-README.md','reports/study/completed-compressive160/README.md')
    fresh=read(HERE/'verification.json')
    prior=read(STUDY/'completed-compressive160/verification.json')
    assert fresh['status']==prior['status']=='passed'
    assert fresh['saved_fields_checked']==5 and prior['saved_fields_checked']==15
    assert fresh['source_hashes']==prior['source_hashes']
    results={p.parent.name:read(p) for p in (STUDY/'runs').glob('*/result.json')}
    published={k:v for k,v in results.items() if v['job']['case'] in ('aligned','compressive') and v['status']=='complete'}
    assert len(published)==32 and all(x in published for x in IDS)
    checks=prior['checks']+fresh['checks']
    assert {x['id'] for x in checks}==set(IDS)
    for check in checks:
        assert sha(STUDY/'runs'/check['id']/'result.json')==check['result_sha256']
    for name,value in fresh['source_hashes'].items():
        assert sha(STUDY/name)==value
        assert blob((STUDY/name).read_bytes())==base['reports/study/'+name]['sha']
    for rid in published:
        if rid!='compressive-fd4-n160-half':
            assert blob((STUDY/'runs'/rid/'result.json').read_bytes())==base[f'reports/study/runs/{rid}/result.json']['sha']
    q=DEST/'runs/compressive-fd4-n160-half/result.json';q.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(STUDY/'runs/compressive-fd4-n160-half/result.json',q)
    comparisons={k:v for k,v in read(STUDY/'comparisons.json').items() if v['a'] in published and v['b'] in published}
    pairs={}
    for method in ('fourier','fd4'):
        pair=comparisons[f'compressive-{method}-n160-base__compressive-{method}-n160-half']
        assert pair['diagnostics']['Wmax']['last_time']==.4 and len(pair['fields'])==4
        pairs[method]=pair
    now=datetime.now(timezone.utc).isoformat()
    local_complete=sorted(k for k,v in results.items() if v['status']=='complete')
    local_active=sorted(k for k,v in results.items() if v['status']=='running')
    endpoints={rid:published[rid]['rows'][-1] for rid in IDS}
    fc=endpoints[IDS[1]];fd=endpoints[IDS[3]]
    method_gap=abs(fc['Wmax']-fd['Wmax'])/fc['Wmax']
    stats={'snapshot_utc':now,'local_saved_complete':len(local_complete),'published_verified_complete':32,'planned_runs':48,
           'new_completed_run':'compressive-fd4-n160-half','endpoints':endpoints,'timestep_comparisons':pairs,
           'endpoint_W_method_difference_relative_to_fourier_half':method_gap,
           'first_spectral_warning':{r:next((v['t'] for v in published[r]['rows'] if v['high_band_enstrophy_fraction']>.01),None) for r in IDS},
           'first_global_width_warning':{r:next((v['t'] for v in published[r]['rows'] if v['global_width']['minimum_chord_cells']<6),None) for r in IDS},
           'local_complete_awaiting_next_verified_batch':sorted(set(local_complete)-set(published)),
           'local_active':local_active}
    dump(BATCH/'summary.json',stats)
    dump(DEST/'comparisons.json',comparisons)
    shutil.copyfile(HERE/'verification.json',BATCH/'verification.json')
    combined={'status':'passed','checked_utc':now,'newly_remeasured_fields':5,'previously_verified_fields_reused':15,
              'total_retained_fields_covered':20,'source_hashes':fresh['source_hashes'],'checks':checks,
              'prior_verification':'../completed-compressive160/verification.json',
              'scope':'Five new FD4 half-step fields remeasured. Prior 15-field audit reused after checking unchanged result and source hashes. No trajectories rerun.'}
    dump(BATCH/'combined-verification.json',combined)
    for name in ('verify.py','package.py'):
        ast.parse((HERE/name).read_text());shutil.copyfile(HERE/name,BATCH/name)
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    styles=[('#247c91','-','Fourier dt=0.0004'),('#c76825','--','Fourier dt=0.0002'),
            ('#6f58a5',':','FD4 dt=0.0004'),('#218661','-.','FD4 dt=0.0002')]
    for rid,(color,style,label) in zip(IDS,styles):
        rr=published[rid]['rows'];t=[v['t'] for v in rr]
        axes[0,0].plot(t,[v['Wmax'] for v in rr],style,color=color,label=label)
        axes[0,1].plot(t,[100*v['high_band_enstrophy_fraction'] for v in rr],style,color=color)
        axes[0,2].plot(t,[v['global_width']['minimum_chord_cells'] for v in rr],style,color=color)
        axes[1,1].semilogy(t,np.maximum([abs(v['energy_balance_relative']) for v in rr],1e-18),style,color=color)
        axes[1,2].semilogy(t,np.maximum([abs(v['enstrophy_balance_relative']) for v in rr],1e-18),style,color=color)
    for method,color in [('fourier','#247c91'),('fd4','#218661')]:
        a=published[f'compressive-{method}-n160-base']['rows'];b=published[f'compressive-{method}-n160-half']['rows']
        av=np.array([v['Wmax'] for v in a]);bv=np.array([v['Wmax'] for v in b])
        axes[1,0].semilogy([v['t'] for v in b],np.maximum(abs(av-bv)/np.max(abs(bv)),1e-12),color=color,label=method.upper())
    axes[0,1].axhline(1,color='#a23b36',ls=':',label='1% screen')
    axes[0,2].axhline(6,color='#a23b36',ls=':',label='Six-cell screen')
    titles=('Four completed 160-grid controls','Highest retained band','Width at global maximum','Peak-vorticity timestep differences','Absolute energy-budget residual','Absolute enstrophy-budget residual')
    labels=('Maximum vorticity W','Percent of physical enstrophy','Grid cells','Difference / half-step curve maximum','Fraction of initial energy','Fraction of initial native enstrophy')
    for ax,title,label in zip(axes.flat,titles,labels):ax.set(title=title,xlabel='Model time',ylabel=label,xlim=(0,.4));ax.grid(alpha=.18)
    for ax in (axes[0,0],axes[0,1],axes[0,2],axes[1,0]):ax.legend(fontsize=8)
    fig.suptitle('Compressive surroundings: all 160-grid controls through 0.40',fontsize=15)
    fig.savefig(BATCH/'control-checks.png',dpi=150);fig.savefig(BATCH/'control-checks.svg');plt.close(fig)
    f_err=100*pairs['fourier']['diagnostics']['Wmax']['max_curve_difference_relative_to_reference_peak']
    d_err=100*pairs['fd4']['diagnostics']['Wmax']['max_curve_difference_relative_to_reference_peak']
    d_vel=100*pairs['fd4']['fields'][-1]['relative_velocity_L2_common_modes']
    table=['| Method | Step | Final W | Final I | High-band enstrophy | Global width (cells) | Central width (cells) |',
           '| --- | --- | --- | --- | --- | --- | --- |']
    for rid in IDS:
        z=endpoints[rid];j=published[rid]['job']
        table.append(f"| {j['method']} | {j['dt']:.4f} | {z['Wmax']:.6f} | {z['I']:.6f} | {100*z['high_band_enstrophy_fraction']:.4f}% | {z['global_width']['minimum_chord_cells']:.4f} | {z['core_width']['minimum_chord_cells']:.4f} |")
    report=f'''# Compressive vortex: all four 160-grid controls complete

All four controls reached model time **0.40**. This completes all **16 compressive-case runs** and brings the verified published matrix to **32/48** (16 aligned and 16 compressive). The separate departure case is progressing locally.

![Completed controls and budgets](control-checks.png)

## Measured comparisons

- Halving the timestep changes the peak-vorticity curve by **{f_err:.6f}% for Fourier** and **{d_err:.6f}% for FD4**. Each percentage is the maximum absolute curve difference divided by its half-step curve maximum.
- The FD4 endpoint velocity L2 difference is **{d_vel:.6f}%** over the common retained modes.
- Both Fourier timesteps exceed the 1% high-band-enstrophy screen from **t=0.28**; both FD4 timesteps exceed it from **t=0.30**. Endpoint fractions are **{100*fc['high_band_enstrophy_fraction']:.4f}%** and **{100*fd['high_band_enstrophy_fraction']:.4f}%** respectively for the half-step runs.
- At the endpoint, the half-step methods differ in W by **{100*method_gap:.4f}%**, relative to Fourier. Their global maxima occur at different locations. Fourier's global width is **{fc['global_width']['minimum_chord_cells']:.4f} cells** and FD4's is **{fd['global_width']['minimum_chord_cells']:.4f} cells**. Their central-region widths are **{fc['core_width']['minimum_chord_cells']:.4f}** and **{fd['core_width']['minimum_chord_cells']:.4f} cells**.
- The FD4 half-step absolute relative energy-budget residual is **{abs(fd['energy_balance_relative']):.6g}**; its native-enstrophy-budget residual is **{abs(fd['enstrophy_balance_relative']):.6g}**.

{chr(10).join(table)}

Timestep agreement is strong, while the spatial-resolution screens fail. The difference in global peak locations matters when interpreting the width and W comparisons. These calculations do not establish spatial convergence, singularity formation, or central-core collapse. All recorded warnings remain in the result files.

## Verification

The new FD4 half-step result was checked against **five retained fields** at t=0, 0.1, 0.2, 0.3 and 0.4. Independent physical-space sums reproduce maximum vorticity, energy, native enstrophy and Fourier-curl enstrophy. Frozen diagnostics reproduce endpoint widths, strain and spectra. Source hashes, settings, finite arrays, zero mean, saved accumulators and budget gates were checked.

The earlier audit of the other **15 retained fields** is reused after confirming unchanged source and result hashes. The combined record covers **20 fields across four runs**. Recorded stage-integrated I and production/loss accumulations were checked against saved checkpoints; the trajectories were not rerun.

[New five-field verification](verification.json) · [Combined verification](combined-verification.json) · [Measurements](summary.json) · [Completed comparisons](../comparisons.json) · [Earlier three-run report](../completed-compressive160/README.md)

## Preserved model and scope

This separate study uses its prescribed smooth periodic vortex start, cube side 6, viscosity 0.001, zero force and SSP RK3. FD4 discretizes transport and viscosity separately while sharing FFT pressure/filter infrastructure and the reported Fourier-curl W diagnostic. The original matched-stretch study retains its supplied field, nonsmooth periodic joins and initial maxima outside the central tube.

[Solver](../solver.py) · [Initial construction](../initial_design.py) · [Protocol](../protocol.json) · [Runner](../run_study.py) · [Study index](../README.md)

## Reproduction

With Python 3.12, NumPy and SciPy, run the existing study runner from its directory:

    python run_study.py --case compressive --n 160 --method fd4 --half

The runner protects active work with a process lock. Full restart fields stay in the local workspace because of their size; SHA-256 hashes are retained in the verification records. With those fields available, the new audit is:

    python completed-compressive160-final/verify.py

The package builder uses the recorded remote-base files and local results; it performs no numerical evolution.
'''
    write(BATCH/'README.md',report)
    update='> Update: all four 160-grid controls are now complete. See the [verified four-run report](../completed-compressive160-final/README.md). The three-run report below preserves its earlier publication snapshot.\n\n'
    write(DEST/'completed-compressive160/README.md',update+old_batch)
    rows=['| Run | Status | Time | Largest saved W | Last W | Last I |','| --- | --- | --- | --- | --- | --- |']
    for rid,r in sorted(published.items()):
        z=r['rows'][-1];rows.append(f"| {rid} | complete | 0.40 | {max(v['Wmax'] for v in r['rows']):.5f} | {z['Wmax']:.5f} | {z['I']:.5f} |")
    footer=old_study[old_study.index('These are different starting fields'):].strip()
    footer=footer.replace('[New three-run verification]','[Earlier three-run verification]')
    main=f'''# Vortex comparison: current numerical results

Snapshot: {now}. **{len(local_complete)}/48 runs have completed locally; 32 verified completed records are published below.** All 16 aligned and all 16 compressive runs are complete and published. The departure case is progressing; its new results await their verified batch.

## All four 160-grid compressive controls

[Completed controls, verification and resolution findings](completed-compressive160-final/README.md). Timestep differences are small for both methods; both still fail the spectral screen. The Fourier and FD4 global maxima differ in location and measured width.

![Completed controls](completed-compressive160-final/control-checks.png)

{chr(10).join(rows)}

{footer}

[New FD4 verification and combined four-run audit](completed-compressive160-final/combined-verification.json).
'''
    write(DEST/'README.md',main)
    old_line=next(x for x in old_index.splitlines() if '[Three vortex surroundings]' in x)
    new_line=f'- **[Three vortex surroundings](study/README.md): {len(local_complete)}/48 completed locally; 32 verified results published.** All aligned and compressive runs are published, including [four 160-grid controls](study/completed-compressive160-final/README.md). Departure runs are progressing; the older interactive display remains dated separately below.'
    write(PUB/'reports/README.md',old_index.replace(old_line,new_line,1))
    dump(DEST/'report-status.json',{'completed':len(local_complete),'published_verified_complete':32,'planned':48,'study_status':'running','snapshot_utc':now,'counts_derived_from':'local result.json completion status; published table contains verified aligned and compressive records','comparisons':len(comparisons),'new_batch':'completed-compressive160-final/README.md'})
    dump(DEST/'task-state.json',{'status':'running','snapshot_utc':now,'planned':48,'completed':[{'id':k,'status':'complete'} for k in local_complete],'published_verified_complete':sorted(published),'active':local_active,'queued':[j['id'] for j in read(STUDY/'protocol.json')['jobs'] if j['id'] not in results],'source_hashes':fresh['source_hashes'],'scope':'Current local status snapshot; departure results await their verified publication batch. Process IDs and checkpoints remain local.'})
    manifest=[]
    for f in sorted(PUB.rglob('*')):
        if f.is_file():
            b=f.read_bytes();manifest.append({'path':f.relative_to(PUB).as_posix(),'sha256':hashlib.sha256(b).hexdigest(),'git_blob_sha':blob(b),'bytes':len(b)})
    dump(HERE/'publication-payload.json',manifest)
    print(json.dumps({'files':len(manifest),'published_complete':32,'local_complete':len(local_complete),'new_fields_verified':5,'FD4_W_curve_difference_percent':d_err,'method_endpoint_W_difference_percent':100*method_gap}))
if __name__=='__main__':main()

