"""Verify saved data and physical invariants; failed research screens are expected."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.signal import resample
from normal_modes import read_hessian,analyze_modes
from run_fluid import properties,exact_startup
from vendor.ns_solver import Solver

def read(path):return json.loads(path.read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(root):
    tests=[]
    def check(condition,name):
        if not condition:raise AssertionError(name)
        tests.append(name)
    provenance=read(root/'provenance.json')
    # A rerun folder contains results; the unchanged source remains beside this script.
    vendor=Path(__file__).resolve().parent/'vendor/ns_solver.py'
    check(digest(vendor)==provenance['reused_file_sha256'],'Reused Navier-Stokes solver matches the pinned repository file')
    records=read(root/'vibration-final/runs.json');cfg=read(root/'protocol.json')
    check(len(records)==12,'All 12 primary vibration executions completed')
    for r in records:
        folder=root/'vibration-final'/r['tag'];lines=(folder/'xtbopt.xyz').read_text().splitlines()[2:]
        xyz=np.array([[float(v) for v in line.split()[1:4]] for line in lines if line.strip()])
        H=read_hessian(folder/'hessian');mH=cfg['mass_u']['H1' if r['isotope']=='H2O' else 'D2'];masses=[cfg['mass_u']['O16'],mH,mH]
        frequencies,allm,_=analyze_modes(H,xyz,masses)
        check(np.allclose(frequencies,r['frequency_cm1'],rtol=0,atol=1e-8),'Recomputed Hessian modes: '+r['tag'])
        check(np.all(frequencies>0),'Positive internal modes: '+r['tag'])
    md=read(root/'molecular-data/summary.json')
    check(md['complete_main_seeds']==[1,2,3],'All independent main seeds completed')
    check(len(md['runs'])==24,'All 24 planned molecular cases completed')
    for r in md['runs']:
        folder=root/'molecular-data'/r['tag'];trajectory=np.load(folder/'trajectory.npz');x=trajectory['positions_nm']
        check(x.shape==(round(r['production_ps']/r['sample_ps']),r['beads'],r['n'],3,3) and np.isfinite(x).all(),'Finite molecular trajectory and expected shape: '+r['tag'])
        if r['mobility_NVE_ps']:
            check((folder/'mobility.npz').is_file(),'Mobility branch saved: '+r['tag'])
    pp=read(root/'fluid-data/properties.json')
    for T,v in pp.items():
        p=properties(float(T))
        check(all(abs(p[i][k]/v[i][k]-1)<1e-12 for i in p for k in p[i]),'Recomputed constitutive properties: '+T)
    channels=read(root/'fluid-data/channel_results.json')
    for r in channels:
        data=np.load(root/'fluid-data'/f"channel-T{r['T_K']}-{r['material']}-G{r['G_Pam']}.npz")
        x=data['y_m'];u=data['profiles_ms'][-1]
        exact=exact_startup(x,1.,r['rho'],r['mu'],r['G_Pam'],r['H_m'])
        check(np.linalg.norm(u-exact)/np.linalg.norm(exact)<1e-3,'Channel field matches analytic solution: '+str((r['T_K'],r['material'],r['G_Pam'])))
        check(r['energy_relative_residual']<1e-3,'Channel energy conservation: '+str((r['T_K'],r['material'],r['G_Pam'])))
    check(read(root/'fluid-data/exact_vortex_check.json')['relative_error']<1e-9,'Exact 2D vortex independent solver verification')
    budget=read(root/'channel-error-budget.json')['comparisons']
    check(len(budget)==6 and max(r['reconstructed_vs_recorded_error_difference'] for r in budget)<1e-11,'Independent discrete-eigenmode reconstruction of channel errors')
    check(len(read(root/'fluid-stress/results.json'))==12,'All stronger-flow stress cases completed')
    check(len(read(root/'fluid-refinement/results.json'))==6,'All bounded finer-grid and half-step cases completed')
    refinement=read(root/'fluid-refinement/convergence.json')
    for iso in ['H2O','D2O']:
        coarse=np.load(root/'fluid-stress'/f'{iso}-U0.1-n48.npz')['final_hat']
        fine=np.load(root/'fluid-refinement'/f'{iso}-n64-dt0.01.npz')['final_hat']
        u=Solver(48,1,workers=1).real(coarse);v=Solver(64,1,workers=1).real(fine)
        for axis in [1,2,3]:u=resample(u,64,axis=axis)
        direct=np.linalg.norm(u-v)/np.linalg.norm(v)
        recorded=next(r['relative_field_change'] for r in refinement if r['isotope']==iso and r['kind']=='grid' and r['fine_n']==64)
        check(abs(direct-recorded)<1e-12,'Physical-space interpolation independently verifies grid comparison: '+iso)
    check(read(root/'molecular-summary.json')['bulk_diffusion_validated'] is False,'Report retains unvalidated transport status')
    check((root/'report.html').is_file(),'Self-contained report exists')
    return tests

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--write-manifest',action='store_true');a=ap.parse_args();root=a.data
    manifest=root/'checksums.json'
    if manifest.exists() and not a.write_manifest:
        saved=read(manifest)
        for name,expected in saved.items():
            if digest(root/name)!=expected:raise AssertionError('Checksum mismatch: '+name)
        print('Verified',len(saved),'file checksums.')
    tests=verify(root)
    out=dict(status='PASS: execution integrity and implementation checks, not every research hypothesis',checks_passed=len(tests),checks=tests)
    if a.write_manifest:
        (root/'verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
        files={p.relative_to(root).as_posix():digest(p) for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p!=manifest}
        manifest.write_text(json.dumps(files,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
