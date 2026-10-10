"""Independent physical checks, including PIMD harmonic equilibrium variance."""
import argparse,json
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit
from scipy.constants import c, hbar, Avogadro
from normal_modes import analyze_modes,read_hessian

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);ap.add_argument('--data',type=Path,default=Path(__file__).resolve().parent);a=ap.parse_args()
    root=a.data
    records=json.loads((root/'vibration-final/runs.json').read_text());cfg=json.loads((Path(__file__).resolve().parent/'protocol.json').read_text())
    checks={};mass=np.array([cfg['mass_u']['O16'],cfg['mass_u']['H1'],cfg['mass_u']['H1']])
    data=np.load(root/'vibration-final/gfn2-H2O-start0-step0.00125/normal_modes.npz')
    h=data['hessian'];xyz=data['xyz_A'];f,_,_=analyze_modes(h,xyz,mass)
    g,_,_=analyze_modes(h,xyz,4*mass)
    checks['all_masses_times_four_frequency_halving_error_cm1']=float(np.max(abs(g-f/2)))
    d=np.load(root/'vibration-final/gfn2-D2O-start0-step0.00125/normal_modes.npz')
    checks['H_D_electronic_Hessian_max_difference']=float(np.max(abs(h-d['hessian'])))
    checks['six_external_modes_max_absolute_cm1']=max(max(abs(v) for v in r['all_modes_cm1'][:6]) for r in records)
    predictions={i:np.array(next(r['frequency_cm1'] for r in records if r['method']=='GFN2-xTB' and r['isotope']==i and r['hessian_step_bohr']==.00125)) for i in ['H2O','D2O']}
    obs={i:np.array(cfg[i+'_observed_cm-1']) for i in predictions}
    maxerr=max(float(np.max(abs(predictions[i]/obs[i]-1))) for i in obs)
    shifterr=float(np.max(abs(predictions['D2O']/predictions['H2O']-obs['D2O']/obs['H2O'])))
    spread=max(float(np.max(np.ptp([r['frequency_cm1'] for r in records if r['method']=='GFN2-xTB' and r['isotope']==i],axis=0))) for i in obs)
    checks['vibration_gate']=dict(max_absolute_relative_error=maxerr,max_isotope_shift_fraction_error=shifterr,max_numerical_spread_cm1=spread,passed=maxerr<.05 and shifterr<.02 and spread<1.)
    assert checks['vibration_gate']['passed']
    assert checks['all_masses_times_four_frequency_halving_error_cm1']<1e-8
    assert checks['H_D_electronic_Hessian_max_difference']<1e-8
    # 3D isotropic harmonic oscillator has a known finite-bead covariance.
    # This checks OpenMM's temperature convention and the sampling algorithm.
    oscillator=[];T=300.;m=1.00782503223;omega=2*np.pi*c*100*3700*1e-12
    k=m*omega**2;kBT=.008314462618*T;hb=hbar*Avogadro*1e9
    exact=hb/(2*m*omega)/np.tanh(hb*omega/(2*kBT))
    for beads in [1,16,32,64]:
        system=mm.System();system.addParticle(m)
        force=mm.CustomExternalForce('k/2*(x*x+y*y+z*z)');force.addGlobalParameter('k',k);force.addParticle(0,[]);system.addForce(force)
        it=mm.RPMDIntegrator(beads,T,5.,.00025);it.setRandomNumberSeed(2200+beads)
        ctx=mm.Context(system,it,mm.Platform.getPlatformByName('Reference'))
        rng=np.random.default_rng(9200+beads)
        for b in range(beads):it.setPositions(b,np.zeros((1,3)));it.setVelocities(b,rng.normal(size=(1,3))*np.sqrt(beads*kBT/m))
        it.step(40000)
        samples=[]
        for j in range(2000):
            it.step(20)
            samples.append(np.mean([np.mean(np.asarray(it.getState(b,getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer))**2) for b in range(beads)]))
        omegaP=beads*kBT/hb
        finite=kBT*np.sum(1/(k+4*m*omegaP**2*np.sin(np.pi*np.arange(beads)/beads)**2))
        measured=float(np.mean(samples));record=dict(beads=beads,measured_variance_nm2=measured,finite_bead_exact_nm2=finite,continuum_quantum_exact_nm2=exact,classical_exact_nm2=kBT/k,relative_sampling_error=measured/finite-1)
        oscillator.append(record);print(record,flush=True)
        assert abs(measured/finite-1)<.12
        del ctx,it
    checks['harmonic_sampling']=oscillator
    a.out.write_text(json.dumps(checks,indent=2))

if __name__=='__main__':main()
