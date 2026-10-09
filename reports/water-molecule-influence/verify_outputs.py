"""Check saved data, parameter identity, symmetry, and numeric consistency."""
from pathlib import Path
import json
import numpy as np
import openmm as mm
from openmm import unit

D=Path(__file__).resolve().parent/'data'
system=mm.XmlSerializer.deserialize((D/'water_system.xml').read_text())
nonbonded=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
parameters=[]
for i in range(system.getNumParticles()):
    q,sigma,epsilon=nonbonded.getParticleParameters(i)
    parameters.append([system.getParticleMass(i).value_in_unit(unit.dalton),
        q.value_in_unit(unit.elementary_charge),sigma.value_in_unit(unit.nanometer),
        epsilon.value_in_unit(unit.kilojoule_per_mole)])
parameters=np.array(parameters).reshape(216,3,4)
assert np.array_equal(parameters,np.broadcast_to(parameters[0],parameters.shape))
ta=np.load(D/'toy_arrays.npz')
assert np.max(np.abs(ta['response']-ta['response'].T))<1e-12
assert np.max(np.abs(ta['response'].sum(axis=0)-1))<1e-12
ca=np.load(D/'md_control_arrays.npz')
J=ca['base_jacobian']
pair=np.sqrt((J**2).sum(axis=(1,3)));np.fill_diagonal(pair,0)
assert np.allclose(pair.sum(axis=0),ca['base_scores'],rtol=1e-14)
assert np.allclose(ca['permuted_scores'],ca['base_scores'][ca['permutation']],rtol=1e-8)
checks={'identical_mass_charge_LJ_parameters_all_216_waters':True,
        'particle_count':system.getNumParticles(),'constraint_count':system.getNumConstraints(),
        'score_reconstruction_max_absolute_error':float(np.max(abs(pair.sum(axis=0)-ca['base_scores']))),
        'all_saved_arrays_finite':True,'replicas':[]}
for rep in range(3):
    data=np.load(D/f'md_replica_{rep}.npz')
    assert data['positions_nm'].shape==(500,216,3,3)
    assert data['force_scores'].shape==(10,216)
    for key in data.files:assert np.isfinite(data[key]).all()
    xyz=data['positions_nm']
    bonds=np.linalg.norm(xyz[:,:,1:,:]-xyz[:,:,0:1,:],axis=3)
    hh=np.linalg.norm(xyz[:,:,1,:]-xyz[:,:,2,:],axis=2)
    checks['replicas'].append({'replica':rep,
        'max_OH_length_deviation_nm':float(np.max(abs(bonds-.09572))),
        'max_HH_length_deviation_nm':float(np.max(abs(hh-2*.09572*np.sin(104.52*np.pi/360))))})
    assert checks['replicas'][-1]['max_OH_length_deviation_nm']<1e-5
assert system.getNumConstraints()==648
(D/'verification.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
