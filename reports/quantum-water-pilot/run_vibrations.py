"""Execute xTB independently for all specified monomer validation cases.

Usage: python run_vibrations.py --xtb /path/to/xtb --out new-run-directory
The output directory must not already exist. No frequencies are fitted.
"""
import argparse, hashlib, json, os, re, subprocess, time
from pathlib import Path
import numpy as np
from normal_modes import analyze_modes, read_hessian

ROOT = Path(__file__).resolve().parent

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xtb', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    exe = args.xtb.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    cfg = json.loads((ROOT/'protocol.json').read_text())
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    version = subprocess.run([str(exe),'--version'],capture_output=True,text=True,encoding='utf-8',errors='replace',env=env)
    (args.out/'version.txt').write_text(version.stdout+version.stderr,encoding='utf-8')
    records = []
    # Three independent starts at default Hessian displacement, plus two
    # refinements at the first start. GFN1 is a declared secondary diagnostic.
    cases = [(2,k,0.005) for k in range(3)]+[(2,0,h) for h in (0.0025,0.00125)]+[(1,0,0.0025)]
    for method,start,step in cases:
        for isotope in ('H2O','D2O'):
            tag = f'gfn{method}-{isotope}-start{start}-step{step}'
            folder = args.out/tag
            folder.mkdir()
            geom = cfg['starting_configurations'][start]
            r,theta = geom['r_angstrom'], np.deg2rad(geom['angle_deg']/2)
            coords = np.array([[0,0,0],[r*np.sin(theta),0,r*np.cos(theta)],[-r*np.sin(theta),0,r*np.cos(theta)]])
            (folder/'water.xyz').write_text('3\ninitial water in angstrom\n'+''.join(f'{a} {x:.14f} {y:.14f} {z:.14f}\n' for a,(x,y,z) in zip(('O','H','H'),coords)))
            m = cfg['mass_u']['H1' if isotope=='H2O' else 'D2']
            control = f'$hess\n step={step}\n sccacc=0.01\n isotope: 1,{cfg["mass_u"]["O16"]}\n isotope: 2,{m}\n isotope: 3,{m}\n$end\n'
            (folder/'settings.inp').write_text(control)
            cmd = [str(exe),'water.xyz','--gfn',str(method),'--ohess','extreme','--acc','0.01','--input','settings.inp']
            tick = time.perf_counter()
            result = subprocess.run(cmd,cwd=folder,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180)
            log = result.stdout+'\n'+result.stderr
            (folder/'stdout.txt').write_text(log,encoding='utf-8')
            if result.returncode or 'normal termination' not in log:
                raise RuntimeError(f'Failed {tag}; see stdout.txt')
            section = log.split('projected vibrational frequencies')[-1].split('reduced masses')[0]
            freqs = np.array([float(v) for line in section.splitlines() if 'eigval :' in line for v in line.split(':',1)[1].split()])
            if len(freqs)!=9: raise RuntimeError(f'Expected 9 modes for {tag}, got {freqs}')
            optimized = (folder/'xtbopt.xyz').read_text().splitlines()
            xyz = np.array([[float(v) for v in ln.split()[1:4]] for ln in optimized[2:5]])
            bonds = xyz[1:]-xyz[0]
            lengths = np.linalg.norm(bonds,axis=1)
            angle = np.rad2deg(np.arccos(np.dot(*bonds)/np.prod(lengths)))
            hessian = read_hessian(folder/'hessian')
            corrected, all_modes, vectors = analyze_modes(hessian,xyz,[cfg['mass_u']['O16'],m,m])
            # This build resets isotope masses internally: retain its output as
            # a diagnostic, and use only explicitly mass-weighted eigenvalues.
            builtin_check,_,_ = analyze_modes(hessian,xyz,[15.99940492,1.00794075,1.00794075])
            if np.max(abs(np.sort(builtin_check)-freqs[-3:]))>0.15:
                raise RuntimeError('Independent normal modes do not reproduce default xTB masses')
            np.savez(folder/'normal_modes.npz',hessian=hessian,xyz_A=xyz,masses_u=[cfg['mass_u']['O16'],m,m],eigenvectors=vectors)
            rec = dict(tag=tag,method=f'GFN{method}-xTB',isotope=isotope,start=start,hessian_step_bohr=step,frequency_cm1=corrected.tolist(),all_modes_cm1=all_modes.tolist(),builtin_untrusted_frequency_cm1=freqs[-3:].tolist(),builtin_reproduction_max_cm1=float(np.max(abs(np.sort(builtin_check)-freqs[-3:]))),bond_lengths_A=lengths.tolist(),angle_deg=float(angle),seconds=time.perf_counter()-tick,command=cmd)
            records.append(rec)
            print(json.dumps(rec),flush=True)
            (args.out/'runs.json').write_text(json.dumps(records,indent=2))
    (args.out/'execution.json').write_text(json.dumps(dict(executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),run_count=len(records)),indent=2))

if __name__=='__main__': main()
