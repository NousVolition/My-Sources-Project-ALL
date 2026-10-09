"""Analytic checks of the actual OpenMM drag-split integrator."""
from pathlib import Path
import json
import numpy as np
import openmm as mm
from openmm import unit
from scipy.linalg import expm
from damped_water import damped_integrator


def simulate(gamma, dt, spring):
    system = mm.System()
    system.addParticle(1)
    if spring:
        force = mm.CustomExternalForce('0.5*x*x')
        force.addParticle(0, [])
        system.addForce(force)
    integrator = damped_integrator(dt, gamma)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
    context.setPositions([[1., 0, 0]]*unit.nanometer)
    context.setVelocities([[0. if spring else 1., 0, 0]]*unit.nanometer/unit.picosecond)
    values = []
    for _ in range(round(4/dt)):
        integrator.step(1)
        state = context.getState(getPositions=True, getVelocities=True)
        values.append([state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)[0, 0],
                       state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)[0, 0]])
    del context, integrator
    return np.arange(1, len(values)+1)*dt, np.array(values)


def main():
    rows = []
    for gamma in [0., .5, 2., 5.]:
        errors = []
        for dt in [.02, .01]:
            t, values = simulate(gamma, dt, True)
            exact = np.array([expm(np.array([[0, 1], [-1, -gamma]])*s) @ [1, 0] for s in t])
            errors.append(float(np.max(abs(values[:, 0]-exact[:, 0]))))
        assert errors[0]/errors[1] > 3.8
        rows.append({'gamma': gamma, 'max_position_error_dt_0_01': errors[1], 'error_ratio': errors[0]/errors[1]})
    t, values = simulate(5, .01, False)
    exact_v = np.exp(-5*t)
    discrete_x = 1+.01*np.exp(-.025)*np.expm1(-5*t)/np.expm1(-.05)
    x_error, v_error = float(np.max(abs(values[:, 0]-discrete_x))), float(np.max(abs(values[:, 1]-exact_v)))
    assert max(x_error, v_error) < 1e-12
    results = {'spring': rows, 'free_drag_discrete_position_error_nm': x_error,
               'free_drag_velocity_error_nm_ps': v_error, 'passed': True}
    path = Path(__file__).resolve().parent/'data'/'integrator_validation.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
