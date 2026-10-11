"""Convert supplied physical component values; does not infer missing measurements."""
import argparse,json,math
from scipy.constants import elementary_charge as e,hbar,Boltzmann as kB
def parameters(ic1_A,ic2_A,shunt_ohm,load_ohm,temperature_K,capacitance_F=0):
    if min(ic1_A,ic2_A,shunt_ohm,load_ohm)<=0 or temperature_K<0 or capacitance_F<0:
        raise ValueError('Critical currents/resistances must be positive; T and C nonnegative')
    if not all(math.isfinite(x) for x in [ic1_A,ic2_A,shunt_ohm,load_ohm,temperature_K,capacitance_F]):raise ValueError('Finite measured inputs required')
    ic=(ic1_A+ic2_A)/2
    return {'Ic_reference_A':ic,'alpha':shunt_ohm/load_ohm,'critical_current_coefficients':[ic1_A/ic,ic2_A/ic],'Gamma':2*e*kB*temperature_K/(hbar*ic),'beta':2*e*ic*shunt_ohm**2*capacitance_F/hbar,'time_unit_seconds':hbar/(2*e*ic*shunt_ohm),'voltage_unit_V':ic*shunt_ohm,'assumptions':'Equal individual shunt resistances and capacitances; all resistors at the supplied T; classical white thermal noise. Instrument bandwidth, colored noise and quantum effects require separate validation.'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('ic1_A','ic2_A','shunt_ohm','load_ohm','temperature_K'):p.add_argument('--'+name,type=float,required=True)
    p.add_argument('--capacitance_F',type=float,default=0)
    print(json.dumps(parameters(**vars(p.parse_args())),indent=2))
