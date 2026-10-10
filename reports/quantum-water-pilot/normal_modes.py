"""Mass-weight an electronic Cartesian Hessian independently of xTB's isotope I/O."""
import numpy as np
from scipy.constants import physical_constants, c

CONVERSION = np.sqrt(physical_constants['Hartree energy'][0]/physical_constants['Bohr radius'][0]**2/physical_constants['atomic mass constant'][0])/(2*np.pi*c*100)

def analyze_modes(hessian, xyz, masses):
    m = np.repeat(masses,3)
    weighted = hessian/np.sqrt(np.outer(m,m))
    eig, vec = np.linalg.eigh(weighted)
    all_freq = np.sign(eig)*np.sqrt(abs(eig))*CONVERSION
    # The three nonzero internal modes; classify stretches by the signs of
    # their two O-H bond-length changes rather than blindly sorting modes.
    internal = list(np.argsort(eig)[-3:])
    bend = internal.pop(0)
    bonds = xyz[1:]-xyz[0]
    unit = bonds/np.linalg.norm(bonds,axis=1)[:,None]
    assignment = {'bend':bend}
    for k in internal:
        displacement = (vec[:,k]/np.sqrt(m)).reshape(3,3)
        changes = np.sum((displacement[1:]-displacement[0])*unit,axis=1)
        assignment['symmetric stretch' if np.prod(changes)>0 else 'antisymmetric stretch'] = k
    if len(assignment)!=3: raise ValueError('Ambiguous stretch assignment')
    idx = [assignment[name] for name in ('bend','symmetric stretch','antisymmetric stretch')]
    return all_freq[idx],all_freq,vec[:,idx]

def read_hessian(path):
    values = [float(v.replace('D','E')) for line in path.read_text().splitlines() if not line.strip().startswith('$') for v in line.split()]
    h = np.array(values).reshape(9,9)
    if np.max(abs(h-h.T))>1e-7: raise ValueError('Hessian is not symmetric')
    return (h+h.T)/2
