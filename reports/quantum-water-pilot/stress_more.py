"""Additional, explicitly post-primary stress tests; no refitting of the model."""
import argparse,json
from pathlib import Path
import numpy as np
from normal_modes import read_hessian,CONVERSION
from run_fluid import properties,dump
from vendor.ns_solver import Solver,difference

def isotope(root):
    records=json.loads((root/'vibration-final/runs.json').read_text())
    ref=next(r for r in records if r['isotope']=='H2O' and r['method']=='GFN2-xTB' and r['hessian_step_bohr']==.00125)
    H=read_hessian(root/'vibration-final'/ref['tag']/'hessian')
    O=15.99491461957;h=1.00782503223;d=2.01410177812
    def freq(masses):
        m=np.repeat(masses,3)
        return np.sqrt(np.linalg.eigvalsh(H/np.sqrt(np.outer(m,m)))[-3:])*CONVERSION
    predicted=freq([O,h,d]);observed=np.array([1402.20,2726.73,3707.47])
    values=np.array([freq([O,x,x]) for x in np.linspace(h,16,101)])
    result=dict(purpose='Post-primary HDO isotope generalization; same Hessian, no fit. Historical method-training overlap not audited.',source='https://webbook.nist.gov/cgi/cbook.cgi?ID=C14940637&Mask=800',mode_order=['bend','OD stretch','OH stretch'],predicted_cm1=predicted,observed_IR_cm1=observed,error_percent=100*(predicted/observed-1),max_absolute_error_percent=float(np.max(abs(predicted/observed-1))*100),same_5percent_screen_pass=bool(np.max(abs(predicted/observed-1))<.05),fictitious_mass_sweep_u=[h,16],all_modes_monotonically_decrease=bool(np.all(np.diff(values,axis=0)<0)),H_D_exchange_max_error_cm1=float(np.max(abs(freq([O,h,d])-freq([O,d,h])))))
    dump(root/'isotope-stress.json',result);print(result,flush=True)

def fluid(root):
    out=root/'fluid-stress';out.mkdir(exist_ok=False)
    plan=dict(temperature_K=298.15,U_ms=[.05,.1],grids=[24,32,48],dt_dimensionless=.01,end_seconds=.05,acceptance_relative_field_change=.01,interpretation='Bracket breakdown as initial velocity increases 2.5 to 5 times above the original maximum. If finest comparison exceeds 1%, do not claim a converged prediction or turbulence onset.')
    dump(out/'protocol.json',plan);rows=[];checks=[]
    for U in plan['U_ms']:
        for iso,p in properties(298.15).items():
            previous=None
            for n in plan['grids']:
                s=Solver(n,p['nu']/(U*.001),workers=1);h=s.initial('taylor_green');h0=h.copy();e0=.5*s.inner(h,h);loss=0.;cfl=0.
                duration=.05*U/.001;steps=round(duration/.01)
                for j in range(steps):
                    h,dl,c=s.step(h,duration/steps);loss+=dl;cfl=max(cfl,c)
                obs,_,_=s.observe(h,duration,loss,e0)
                row=dict(isotope=iso,U_ms=U,n=n,Re_scale=U*.001/p['nu'],cfl=cfl,rms_velocity_ms=U*np.sqrt(2*obs['energy']),**obs)
                rows.append(row);np.savez_compressed(out/f'{iso}-U{U}-n{n}.npz',final_hat=h,initial_hat=h0)
                if previous is not None:
                    change=difference(previous[0],h,previous[1]);checks.append(dict(isotope=iso,U_ms=U,coarse_n=previous[1].n,fine_n=n,relative_field_change=change,passed=change<.01))
                previous=(h,s);print('STRESS',iso,U,n,obs['energy_residual'],flush=True)
                dump(out/'results.json',rows);dump(out/'convergence.json',checks)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--fluid',action='store_true');a=ap.parse_args()
    isotope(a.data)
    if a.fluid:fluid(a.data)
