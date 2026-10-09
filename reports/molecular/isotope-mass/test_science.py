"""Mechanistic and analysis checks; results saved to validation.json."""
from pathlib import Path
import json, unittest
import numpy as np
import openmm as mm
from openmm import unit
from simulate import system,context,random_positions,masses,response,dump
from analyze import geometry,episode_data,km_rmst,response_metrics,summarize

RESULTS={}
class ScientificTests(unittest.TestCase):
    def test_mass_changes_only_kinetics(self):
        a,_=system(9,-1);b,_=system(9,4)
        self.assertEqual(a.getNumForces(),b.getNumForces())
        for i in range(a.getNumForces()):self.assertEqual(mm.XmlSerializer.serialize(a.getForce(i)),mm.XmlSerializer.serialize(b.getForce(i)))
        self.assertEqual([a.getConstraintParameters(i) for i in range(a.getNumConstraints())],[b.getConstraintParameters(i) for i in range(b.getNumConstraints())])
        changed=np.flatnonzero(masses(a)!=masses(b)).tolist();self.assertEqual(changed,[13,14])
        q=random_positions(9,0,.4);ca,ia=context(a);cb,ib=context(b)
        for c in [ca,cb]:c.setPositions(q)
        fa=ca.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
        fb=cb.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
        err=float(np.max(abs(fa-fb)));self.assertEqual(err,0);RESULTS['isotope_force_difference']=err;RESULTS['changed_atom_indices']=changed

    def test_hbond_geometry(self):
        x=np.array([[[[0.,0,0],[.09572,0,0],[-.024,.0927,0]],[[.28,0,0],[.28,.09572,0],[.1873,-.024,0]]]])
        d,h,c=geometry(x);self.assertTrue(h[0,0,1]);self.assertTrue(h[0,1,0]);self.assertFalse(h[0,0,0])
        x[0,1]+=np.array([1,0,0]);self.assertFalse(geometry(x)[1].any())
        x[0,1]-=np.array([1,0,0]);x[0,0,1]=[-.09572,0,0];x[0,0,2]=[0,-.09572,0]
        x[0,1,1]=[.37572,0,0];x[0,1,2]=[.28,.09572,0]
        self.assertFalse(geometry(x)[1].any())
        RESULTS['known_hbond_geometry']='passed'

    def test_censoring(self):
        a=np.array([1,1,0,1,1,0,0,1,1],bool)[:,None]
        t,e,left=episode_data(a,.1);np.testing.assert_allclose(t,[.2,.1]);self.assertEqual(e.tolist(),[True,False]);self.assertEqual(left,1)
        self.assertAlmostEqual(km_rmst(np.array([.2,.4]),np.array([True,False]),.6),.4)
        RESULTS['censoring_synthetic']='passed'

    def test_permutation_and_free_targets(self):
        s,_=system();q=random_positions(9,22,.4);c,i=context(s,nve=True);c.setPositions(q);mm.LocalEnergyMinimizer.minimize(c,1,1000)
        q=c.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer);v=np.zeros_like(q)
        perm=np.random.default_rng(4).permutation(9);c2,i2=context(s,nve=True);c.setVelocities(v);c2.setPositions(q.reshape(9,3,3)[perm].reshape(-1,3));c2.setVelocities(v)
        i.step(100);i2.step(100)
        qa=c.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        qb=c2.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer).reshape(9,3,3)[np.argsort(perm)].reshape(-1,3)
        err=float(np.max(abs(qa-qb)));self.assertLess(err,1e-8);RESULTS['permuted_nve_0.1ps_position_error_nm']=err
        for fi in reversed(range(s.getNumForces())):
            if isinstance(s.getForce(fi),mm.NonbondedForce):s.removeForce(fi)
        raw,base,ene,lag=response(s,q,v,sources=[0]);der=(raw[:,:,1]-raw[:,:,0])/.1
        err=float(abs(der[:,:,:,1:]).max());self.assertLess(err,1e-10);RESULTS['no_intermolecular_force_offsource_response']=err
        self.assertEqual(float(abs(raw[:,:,1,0]-raw[:,:,0,0]).max()),0.)

    def test_bootstrap_unit(self):
        r=summarize([1,2,3]);self.assertEqual(r['n'],3);self.assertEqual(r['mean'],2);self.assertGreater(r['ci95'][1],r['ci95'][0]);RESULTS['bootstrap_cluster_unit']='passed'

    def test_heavy_identity_permutation(self):
        perm=np.random.default_rng(8).permutation(9);newheavy=int(np.flatnonzero(perm==4)[0]);a,_=system(9,4);b,_=system(9,newheavy)
        ca,ia=context(a,nve=True);cb,ib=context(b,nve=True);q=random_positions(9,24,.4);ca.setPositions(q);mm.LocalEnergyMinimizer.minimize(ca,1,1000)
        q=ca.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        ca.setVelocities(np.zeros_like(q));cb.setPositions(q.reshape(9,3,3)[perm].reshape(-1,3));cb.setVelocities(np.zeros_like(q));ia.step(100);ib.step(100)
        qa=ca.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        qb=cb.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer).reshape(9,3,3)[np.argsort(perm)].reshape(-1,3)
        err=float(abs(qa-qb).max());self.assertLess(err,1e-8);RESULTS['heavy_molecule_relabel_nve_error_nm']=err

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ScientificTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    dump(Path('validation.json'),dict(tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),checks=RESULTS))
    raise SystemExit(not result.wasSuccessful())
