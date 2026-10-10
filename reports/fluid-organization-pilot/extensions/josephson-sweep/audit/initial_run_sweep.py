"""Run all protocol conditions; raw arrays are saved for independent reanalysis."""
from pathlib import Path
import argparse, hashlib, json, platform, time
import numpy as np
from model import integrate, voltage, power_residual

ROOT = Path(__file__).resolve().parent


def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    protocol_bytes = (ROOT/'protocol.json').read_bytes()
    cfg = json.loads(protocol_bytes)
    signature = hashlib.sha256(protocol_bytes + (ROOT/'model.py').read_bytes()
                               + Path(__file__).read_bytes()).hexdigest()
    if (out/'run_metadata.json').exists():
        old = json.loads((out/'run_metadata.json').read_text())
        if old['signature'] != signature:
            raise RuntimeError('Source/protocol changed; choose a new --out directory')
    started = time.time()
    n = cfg['n_independent_initializations']
    initial = np.random.default_rng(cfg['seed']).uniform(-np.pi, np.pi, (n, 2))
    biases = np.array(cfg['biases']); alphas = np.array(cfg['alphas'])
    # Batch three independently initialized modes at the same sweep index.
    a = np.broadcast_to(alphas[None, :, None], (3, len(alphas), n))
    init = np.broadcast_to(initial, (3, len(alphas), n, 2)).copy()
    np.savez_compressed(out/'initial_conditions.npz', phases=initial, biases=biases, alphas=alphas)
    for s in cfg['settings']:
        name = s['name']; target = out/f'{name}.npz'
        if target.exists():
            print(f'{name}: existing signature-matched recording', flush=True); continue
        previous = init.copy(); starts=[]; settled=[]; ends=[]; means=[]; residual=[]
        phase_records=[]; energy_records=[]; ordered_biases=[]
        for j in range(len(biases)):
            ib = np.broadcast_to(np.array([biases[j], biases[j], biases[-1-j]])[:, None, None], a.shape)
            previous[0] = init[0]
            t, p, d = integrate(previous, ib, a, s['settle']+cfg['base_window'], s['dt'], cfg['record_stride'])
            index = round(s['settle']/cfg['record_stride'])
            starts.append(previous.copy()); settled.append(p[index]); ends.append(p[-1])
            means.append(voltage(p[index], p[-1], cfg['base_window']))
            residual.append(power_residual(p, d, ib))
            phase_records.append(p); energy_records.append(d); ordered_biases.append(ib[:, 0, 0])
            previous = p[-1].copy()
            if j % 4 == 0: print(f'{name}: bias index {j+1}/{len(biases)}', flush=True)
        np.savez_compressed(target, initial=np.array(starts), settled=np.array(settled),
                            endpoint=np.array(ends), voltage=np.array(means), residual=np.array(residual),
                            phases=np.array(phase_records), dissipated=np.array(energy_records),
                            time=t, bias_order=np.array(ordered_biases), alphas=alphas,
                            dt=s['dt'], settle=s['settle'], window=cfg['base_window'])
        metadata={'signature':signature,'protocol_sha256':hashlib.sha256(protocol_bytes).hexdigest(),
                  'python':platform.python_version(),'numpy':np.__version__,
                  'elapsed_seconds_this_invocation':time.time()-started}
        (out/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    wtarget=out/'windows.npz'
    if not wtarget.exists():
        with np.load(out/'both.npz') as z:
            p0=z['settled']; order=z['bias_order']
        ib=np.broadcast_to(order[:, :, None, None],p0.shape[:-1])
        aa=np.broadcast_to(alphas[None,None,:,None],p0.shape[:-1])
        s=cfg['settings'][-1]
        print('Independent window replay from stored post-settling states',flush=True)
        t,p,d=integrate(p0,ib,aa,max(cfg['window_checks']),s['dt'],cfg['record_stride'])
        np.savez_compressed(wtarget,time=t,phases=p,dissipated=d,bias_order=order,
                            alphas=alphas,residual=power_residual(p,d,ib),dt=s['dt'])
    print(f'Completed recordings in {time.time()-started:.1f} seconds this invocation',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',default=ROOT/'data')
    main(parser.parse_args().out)
