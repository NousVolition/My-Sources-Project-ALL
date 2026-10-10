"""Deterministic time-step, field-update, work-balance and phase-analysis checks."""
from pathlib import Path
import json,sys
import numpy as np
import openmm as mm
from openmm import unit
from water_battery import make_system,save,HERE,C
sys.path.insert(0,str(HERE.parent/'response-extension'))
from impulse_response import Probe
from analyze_battery import phase_metrics


def run(dt,update,xml,x,v,side):
    assert abs(round(update/dt)*dt-update)<1e-12,'Field update must contain an integer number of MD steps'
    system=make_system(xml,side)
    # A COM-removal operation is an additional energy sink. Disable it only for
    # this deterministic work-balance audit, and remove initial COM momentum.
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i),mm.CMMotionRemover):system.removeForce(i)
    mass=np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())]).reshape(-1,3)
    v=v.copy();v-=np.sum(mass[:,:,None]*v,axis=(0,1))/mass.sum()
    probe=Probe(system,'OpenCL',dt)
    ctx=probe.context;ctx.setPositions(x.reshape(-1,3));ctx.setVelocities(v.reshape(-1,3));ctx.applyVelocityConstraints(1e-8)
    def energy(groups=None):
        state=ctx.getState(getEnergy=True,**({} if groups is None else dict(groups=groups)))
        return (state.getPotentialEnergy()+(state.getKineticEnergy() if groups is None else 0*unit.kilojoule_per_mole)).value_in_unit(unit.kilojoule_per_mole)
    initial=energy();work=0.
    for step in range(round(.2/update)):
        t=(step+.5)*update;before=energy({1});ctx.setParameter('Ex',.5*np.cos(2*np.pi*t));ctx.setParameter('Ey',.5*np.sin(2*np.pi*t));after=energy({1});work+=after-before
        probe.integrator.step(round(update/dt))
    final=energy();state=ctx.getState(getPositions=True,getVelocities=True)
    xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).reshape(-1,3,3)
    vel=np.asarray(state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond)).reshape(-1,3,3)
    result=dict(dt_ps=dt,field_update_ps=update,initial_total_energy=initial,final_total_energy=final,external_work_kJ_mol=work,work_energy_residual_kJ_mol=final-initial-work,
        relative_work_energy_residual=abs(final-initial-work)/abs(initial))
    del probe.context,probe.integrator
    return xyz,vel,result


def main():
    src=HERE.parent/'response-extension'/'data';xml=(src/'probe_system.xml').read_text();d=np.load(src/'response_replica_0.npz')
    side=json.loads((HERE.parent/'data'/'md_protocol.json').read_text())['box_side_nm'];x=d['positions_nm'];v=d['velocities_nm_per_ps']
    runs=[run(dt,update,xml,x,v,side) for dt,update in [(.001,.01),(.0005,.01),(.00025,.01),(.0005,.005),(.0005,.0025)]]
    comparisons=[]
    for i,label in [(1,'half_MD_step'),(2,'quarter_MD_step'),(3,'half_field_update'),(4,'quarter_field_update')]:
        reference=0 if i<3 else 1
        delta=runs[i][0]-runs[reference][0];delta-=side*np.rint(delta/side)
        comparisons.append(dict(change=label,reference_run_index=reference,position_RMS_nm=float(np.sqrt(np.mean(delta*delta))),velocity_RMS_nm_ps=float(np.sqrt(np.mean((runs[i][1]-runs[reference][1])**2)))))
    # Convergence is measured independently from desired physical outcomes.
    err_coarse=np.linalg.norm(runs[0][0]-runs[2][0]);err_fine=np.linalg.norm(runs[1][0]-runs[2][0]);ratio=float(err_coarse/err_fine)
    field_coarse=np.linalg.norm(runs[1][0]-runs[4][0]);field_fine=np.linalg.norm(runs[3][0]-runs[4][0]);field_ratio=float(field_coarse/field_fine)
    diagnostic=dict(runs=[r[2] for r in runs],comparisons=comparisons,MD_step_error_reduction_ratio=ratio,field_update_error_reduction_ratio=field_ratio)
    save(HERE/'data'/'integration_validation_diagnostic.json',diagnostic)
    print(json.dumps(diagnostic),flush=True)
    assert ratio>2 and field_ratio>2
    assert max(r[2]['relative_work_energy_residual'] for r in runs)<.001
    t=np.arange(0,20,.1);z=np.column_stack([.4*np.cos(np.pi*t),.4*np.sin(np.pi*t),np.zeros(len(t))]);slip=phase_metrics(t,z,1.)
    assert not slip['finite_window_tracking'] and abs(slip['relative_phase_slope']+.5)<1e-8
    gap=np.column_stack([.4*np.cos(2*np.pi*t),.4*np.sin(2*np.pi*t),np.zeros(len(t))]);gap[50:100]=0
    assert not phase_metrics(t,gap,1.)['finite_window_tracking']
    result=dict(horizon_ps=.2,thermostat=False,integrator='constrained velocity Verlet/RATTLE',runs=[r[2] for r in runs],comparisons=comparisons,
        MD_step_error_reduction_ratio=ratio,field_update_error_reduction_ratio=field_ratio,phase_slip_fixture_correct=True,phase_gap_fixture_correct=True,passed=True,
        scope='Short deterministic convergence and external-work balance. Stochastic long trajectories require ensemble comparisons, not pointwise identity.')
    save(HERE/'data'/'integration_validation.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
