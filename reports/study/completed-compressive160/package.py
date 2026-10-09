"""Package verified completed observations and current comparisons for publication."""
from pathlib import Path
from datetime import datetime, timezone
import os, json, shutil, hashlib, ast
import numpy as np

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IDS=('compressive-fourier-n160-base','compressive-fd4-n160-base','compressive-fourier-n160-half')
PUB=HERE/'publication'
DEST=PUB/'reports/study'
BATCH=DEST/'completed-compressive160'

def read(p): return json.loads(p.read_text(encoding='utf-8'))
def dump(p,v): p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def write(p,s): p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s,encoding='utf-8')
def blob(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()

def main():
    audit=read(HERE/'verification.json');assert audit['status']=='passed' and audit['saved_fields_checked']==15
    base={v['path']:v for v in read(HERE/'base-tree.json')['tree'] if v['type']=='blob'}
    for name in audit['source_hashes']:
        assert hashlib.sha256((STUDY/name).read_bytes()).hexdigest()==audit['source_hashes'][name]
        assert blob((STUDY/name).read_bytes())==base['reports/study/'+name]['sha']
    results={p.parent.name:read(p) for p in (STUDY/'runs').glob('*/result.json')}
    complete={k:v for k,v in results.items() if v['status']=='complete'}
    assert len(complete)==31
    for rid in IDS:
        p=STUDY/'runs'/rid/'result.json';r=complete[rid]
        assert hashlib.sha256(p.read_bytes()).hexdigest()==next(v for v in audit['checks'] if v['id']==rid)['result_sha256']
        q=DEST/'runs'/rid/'result.json';q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    for rid in set(complete)-set(IDS):
        # Previously published completed records remain byte-for-byte unchanged.
        assert blob((STUDY/'runs'/rid/'result.json').read_bytes())==base[f'reports/study/runs/{rid}/result.json']['sha']
    comparisons={k:v for k,v in read(STUDY/'comparisons.json').items() if v['a'] in complete and v['b'] in complete}
    pair=comparisons['compressive-fourier-n160-base__compressive-fourier-n160-half']
    assert pair['diagnostics']['Wmax']['last_time']==.4 and len(pair['fields'])==4
    a=complete[IDS[0]];b=complete[IDS[2]];fd=complete[IDS[1]]
    x=np.array([r['t'] for r in b['rows']])
    stats={'completed_runs':31,'planned_runs':48,'newly_published_completed_runs':list(IDS),
           'timestep_comparison':pair,'fourier_endpoint':b['rows'][-1],
           'ordinary_fd4_endpoint':fd['rows'][-1],
           'first_spectral_warning':.28,'first_global_width_warning':.36,
           'width_transition':[r for r in b['rows'] if .339<=r['t']<=.361]}
    dump(BATCH/'summary.json',stats)
    dump(DEST/'comparisons.json',comparisons)
    shutil.copyfile(HERE/'verification.json',BATCH/'verification.json')
    for name in ('verify.py','package.py'):
        ast.parse((HERE/name).read_text(encoding='utf-8'));shutil.copyfile(HERE/name,BATCH/name)
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for rid,label,color,style in ((IDS[0],'Fourier dt=0.0004','#247c91','-'),(IDS[2],'Fourier dt=0.0002','#c76825','--'),(IDS[1],'FD4 dt=0.0004','#777777',':')):
        rows=complete[rid]['rows'];t=[r['t'] for r in rows]
        axes[0,0].plot(t,[r['Wmax'] for r in rows],style,color=color,label=label)
        axes[0,1].plot(t,[100*r['high_band_enstrophy_fraction'] for r in rows],style,color=color)
        axes[0,2].plot(t,[r['global_width']['minimum_chord_cells'] for r in rows],style,color=color)
        axes[1,1].semilogy(t,np.maximum([abs(r['energy_balance_relative']) for r in rows],1e-18),style,color=color)
    for key,label in (('Wmax','Peak vorticity W'),('I','Accumulated I')):
        av=np.array([r[key] for r in a['rows']]);bv=np.array([r[key] for r in b['rows']])
        axes[1,0].semilogy(x,np.maximum(abs(av-bv)/np.max(abs(bv)),1e-12),label=label)
    axes[1,2].plot(x,[r['global_width']['minimum_chord_cells'] for r in b['rows']],label='Global maximum')
    axes[1,2].plot(x,[r['core_width']['minimum_chord_cells'] for r in b['rows']],label='Central region')
    axes[0,1].axhline(1,color='#a23b36',ls=':',label='1% screen')
    for ax in (axes[0,2],axes[1,2]):ax.axhline(6,color='#a23b36',ls=':',label='Six-cell screen')
    titles=('Completed 160-grid trajectories','Highest retained band','Width at global maximum','Fourier timestep differences','Absolute energy balance residual','Half-step widths at two locations')
    labels=('Maximum vorticity','Percent of physical enstrophy','Grid cells','Difference / half-step curve maximum','Fraction of initial energy','Grid cells')
    for ax,title,label in zip(axes.flat,titles,labels):ax.set(title=title,xlabel='Model time',ylabel=label,xlim=(0,.4));ax.grid(alpha=.18)
    for ax in (axes[0,0],axes[0,1],axes[1,0],axes[1,2]):ax.legend(fontsize=8)
    fig.suptitle('Compressive surroundings: completed controls through 0.40',fontsize=15)
    fig.savefig(BATCH/'control-checks.png',dpi=150);fig.savefig(BATCH/'control-checks.svg');plt.close(fig)
    werr=100*pair['diagnostics']['Wmax']['max_curve_difference_relative_to_reference_peak']
    verr=100*pair['fields'][-1]['relative_velocity_L2_common_modes']
    last=b['rows'][-1];transition=stats['width_transition'][-1]
    report=f'''# Compressive vortex: completed 160-grid controls

Three completed runs through model time **0.40**: Fourier at dt=0.0004 and dt=0.0002, and FD4 at dt=0.0004. The FD4 half-step run is still in progress. These records bring the separate vortex matrix to **31/48 completed runs**.

![Control comparisons](control-checks.png)

## Measured results

- Fourier dt-halving changes the maximum-vorticity curve by **{werr:.6f}%**, measured as the maximum absolute difference divided by the half-step curve maximum. Endpoint velocity L2 difference is **{verr:.6f}%** over the common retained modes.
- The Fourier half-step endpoint is W={last['Wmax']:.8f}, I={last['I']:.8f}. Its absolute relative energy-budget residual is {abs(last['energy_balance_relative']):.6g}.
- Both Fourier timesteps first exceed the **1% spectral screen at t=0.28**. At t=0.40, **{100*last['high_band_enstrophy_fraction']:.4f}%** of physical enstrophy lies in the highest retained band.
- Both first fail the global width screen at **t=0.36**. The sampled global maximum moves to a feature **{transition['global_width']['minimum_chord_cells']:.5f} cells** wide; the central-region width remains **{transition['core_width']['minimum_chord_cells']:.5f} cells**. The changing global maximum does not track one material core.

Timestep agreement is strong over this interval, while the spatial-resolution screens fail. These results do not establish spatial convergence, a singularity, or central-core collapse. The warnings are retained in every affected result row.

## Verification and scope

All **15 retained fields** from the three runs were checked at t=0, 0.1, 0.2, 0.3, 0.4. Physical-space sums independently reproduce W, energy and both enstrophies. Endpoint widths, strain and spectra were remeasured with the frozen diagnostic code. Source hashes, zero mean, finite arrays, source/settings provenance and energy-budget gates were checked. Integrated I and budget rates retain the saved RK-stage accumulations; no trajectory was rerun for this audit.

The calculation uses the prescribed smooth periodic vortex construction, cube side 6, viscosity 0.001, zero external force and SSP RK3. This is a separate starting field from matched-stretch. The matched-stretch study retains its original nonsmooth periodic joins and outside-tube initial maxima. FD4 and Fourier share the FFT pressure/filter infrastructure and the reported Fourier-curl W diagnostic.

[Verification](verification.json) · [Detailed measurements](summary.json) · [All completed comparisons](../comparisons.json) · [Study protocol](../protocol.json) · [Study index](../README.md)

## Reproduction

The unchanged [runner](../run_study.py), [solver](../solver.py) and [initial construction](../initial_design.py) reproduce each case using NumPy, SciPy and Python 3.12. Run from the study folder: `python run_study.py --case compressive --n 160 --method fourier` and append `--half` for its timestep control; choose `--method fd4` for the ordinary FD4 run. The runner protects active work with a process lock.

Full restart arrays remain in the local study workspace because of their size; their SHA-256 hashes are in the verification record. With those arrays present, run `python completed-compressive160/verify.py` to repeat the saved-field audit. The original local `package.py` builds this publication using its recorded remote-base snapshot.
'''
    write(BATCH/'README.md',report)
    old=(HERE/'base-study-README.md').read_text(encoding='utf-8')
    rows=['| Run | Status | Time | Largest saved W | Last W | Last I |','| --- | --- | --- | --- | --- | --- |']
    for name,r in sorted(complete.items()):
        lastrow=r['rows'][-1];peak=max(v['Wmax'] for v in r['rows'])
        rows.append(f"| {name} | complete | 0.40 | {peak:.5f} | {lastrow['Wmax']:.5f} | {lastrow['I']:.5f} |")
    now=datetime.now(timezone.utc).isoformat()
    footer=old[old.index('These are different starting fields'):]
    main='# Vortex comparison: current numerical results\n\n'+f'Snapshot: {now}. **31 of 48 runs complete.** All 16 aligned runs and 15 compressive runs are complete. The compressive FD4 half-step run remains active; all 16 departure runs remain queued.\n\n'
    main+='## Newly completed 160-grid controls\n\n[Three completed runs, numerical verification and resolution findings](completed-compressive160/README.md). Both Fourier timesteps fail the spectral screen from 0.28 and global-width screen from 0.36; timestep agreement does not establish spatial convergence.\n\n![Completed controls](completed-compressive160/control-checks.png)\n\n'
    main+='\n'.join(rows)+'\n\n'+footer
    main=main.replace('[Snapshot verification](../numerical-progress-verification.json)','[Earlier snapshot verification](../numerical-progress-verification.json) · [New three-run verification](completed-compressive160/verification.json)')
    write(DEST/'README.md',main)
    index=(HERE/'base-reports-README.md').read_text(encoding='utf-8')
    assert ': 28/48 runs complete.' in index
    index=index.replace('## Earlier numerical snapshot â€” October 8','## Numerical reports and earlier snapshots',1)
    index=index.replace(': 28/48 runs complete.** Updated measurements and comparisons;',': 31/48 runs complete.** Updated October 9 with three verified 160-grid controls and resolution warnings;',1)
    index=index.replace('The counts above describe the saved October 8 snapshot.','The matched-stretch count above describes the saved October 8 snapshot; the vortex count was updated October 9.',1)
    write(PUB/'reports/README.md',index)
    dump(DEST/'report-status.json',{'completed':31,'planned':48,'study_status':'running','snapshot_utc':now,'counts_derived_from':'complete result.json records','comparisons':len(comparisons),'report':'README.md','new_batch':'completed-compressive160/README.md'})
    dump(DEST/'task-state.json',{'status':'running','snapshot_utc':now,'planned':48,'completed':[{'id':k,'status':'complete'} for k in sorted(complete)],'active':['compressive-fd4-n160-half'],'queued':[j['id'] for j in read(STUDY/'protocol.json')['jobs'] if j['id'] not in complete and j['id']!='compressive-fd4-n160-half'],'source_hashes':audit['source_hashes'],'scope':'Publication snapshot derived from saved result status; process identifiers and checkpoints remain local.'})
    manifest=[]
    for path in sorted(PUB.rglob('*')):
        if path.is_file():
            bts=path.read_bytes();manifest.append({'path':path.relative_to(PUB).as_posix(),'sha256':hashlib.sha256(bts).hexdigest(),'git_blob_sha':blob(bts),'bytes':len(bts)})
    dump(HERE/'publication-payload.json',manifest)
    print(json.dumps({'files':len(manifest),'completed_runs':31,'new_completed_runs':3,'W_curve_difference_percent':werr,'endpoint_velocity_L2_percent':verr,'recorded_fields_verified':15}))

if __name__=='__main__':main()
