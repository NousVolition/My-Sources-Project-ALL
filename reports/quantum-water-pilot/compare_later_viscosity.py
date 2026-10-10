import argparse,json
from pathlib import Path
import numpy as np
from run_fluid import properties,dump

def main(root):
    cfg=json.loads((root/'later-viscosity-reference.json').read_text());rows=[]
    for r in cfg['observations']:
        T=r['T_K'];x=r['x_D2O'];p=properties(T)
        def mix(temp):
            q=properties(temp);return 1000*((1-x)*q['H2O']['mu']+x*q['D2O']['mu'])
        observed=float(np.mean(r['readings_mPas']));pred=mix(T)
        dT=(mix(T+.01)-mix(T-.01))/.02
        dx=1000*(p['D2O']['mu']-p['H2O']['mu'])
        sigma=float(np.sqrt((.015*observed)**2+(.15*dT)**2+((.04 if x else 0)*dx)**2))
        rows.append(dict(**r,observed_mean_mPas=observed,predicted_mPas=pred,relative_error_percent=100*(pred/observed-1),approx_uncertainty_mPas=sigma,residual_in_approx_uncertainty_units=(pred-observed)/sigma))
    dump(root/'later-viscosity-comparison.json',dict(protocol=cfg,comparisons=rows,max_absolute_error_percent=max(abs(r['relative_error_percent']) for r in rows)))
    print(rows,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);a=ap.parse_args();main(a.data)
