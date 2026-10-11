import unittest
import numpy as np
from experiment import DrivenSolver,basis,pulse_coefficients,matched_translation
from reference_solver import random_field,interventions
from independent_solver import IndependentSolver

class Controls(unittest.TestCase):
    def test_work_and_strength(self):
        c=np.random.default_rng(3).normal(size=(500,3))*.4
        c=np.concatenate([c,np.zeros((1,3)),np.array([[0,0,1],[0,0,-1]])])
        q,a=pulse_coefficients(c)
        self.assertLess(a['equal_work_error'],1e-13);self.assertLess(a['equal_rms_error'],1e-15)
        self.assertLess(a['orthogonal_linear_work_error'],1e-14)
        self.assertTrue(a['fallback'][-1]);self.assertTrue(a['fallback'][-2])
    def test_basis_orthogonality_and_local_work(self):
        s=DrivenSolver(24,.12,backend='cpu');q=basis(s,[.3,.4,.9]);gram=np.array([[s.inner(a,b) for b in q] for a in q])
        np.testing.assert_allclose(gram,np.eye(3),atol=1e-14)
        h=random_field(s,8)[None];c=np.stack([s.inner(h,q[j]) for j in range(3)],axis=1)
        coefficients,a=pulse_coefficients(c)
        for j in [1,2]:
            impulse=np.einsum('i,icxyz->cxyz',coefficients[0,j],q)
            increment=.5*s.inner(h+impulse,h+impulse)-.5*s.inner(h,h)
            np.testing.assert_allclose(increment,.0002,atol=1e-14)
    def test_forcing_injection(self):
        s=DrivenSolver(24,.04,.2,'cpu');h=random_field(s,2)[None];r,b,_=s.rhs(h,np.zeros_like(h))
        self.assertAlmostEqual(float(b[0]+s.nu*s.inner(h,s.k2*h)[0]),.2,places=14)
    def test_matched_displacement(self):
        s=DrivenSolver(24,.12,backend='cpu');h=random_field(s,5);q=basis(s,[1,2,3],.45,4,True)[2]
        v,checks=interventions(s,np.stack([h+.2*q,h-.2*q]),9,[.13,-.09,.07],4)
        shifted,m=matched_translation(s,v['retained'],np.stack([v['scramble'+str(j)] for j in range(4)]))
        self.assertLess(m['matching_error'],1e-12)
        np.testing.assert_allclose(abs(shifted),abs(v['retained']),atol=1e-12)
    def test_independent_rhs_and_observation(self):
        s=DrivenSolver(24,.04,.2,'cpu');h=np.stack([random_field(s,3),random_field(s,4)])
        independent=IndependentSolver(s.modes_cpu,.04,.2);full=independent.expand(s.pack(h))
        other,rate=independent.rhs(full);reference,rr,_=s.rhs(h,np.zeros_like(h))
        np.testing.assert_allclose(independent.contract(other),s.pack(reference),atol=2e-15)
        np.testing.assert_allclose(rate,rr,atol=2e-15)
        x=np.random.default_rng(2).uniform(0,6.28,(2,4,3));v,j=s.sample(h,x,4);vv,jj=independent.sample(full,x)
        np.testing.assert_allclose(vv,v,atol=3e-14);np.testing.assert_allclose(jj,j,atol=3e-14)
    def test_independent_analytic_shear(self):
        s=DrivenSolver(16,.2,backend='cpu');a=np.arange(16)*s.dx;xx,yy,zz=np.meshgrid(a,a,a,indexing='ij')
        u=np.zeros((1,3,16,16,16));u[0,0]=np.sin(yy);coeff=s.pack(s.hat(u));x=np.array([[[.3,.7,.5]]])
        r=IndependentSolver(s.modes_cpu,.2).integrate(coeff,x,.4)
        factor=-np.expm1(-.08)/.2;true=x.copy();true[...,0]+=np.sin(x[...,1])*factor
        f=np.tile(np.eye(3),(1,1,1,1));f[...,0,1]=np.cos(x[...,1])*factor
        np.testing.assert_allclose(r['positions'],true,atol=2e-12);np.testing.assert_allclose(r['tangent'],f,atol=2e-12)
        self.assertLess(r['energy_budget'],1e-10);self.assertLess(r['volume'],1e-11)

if __name__=='__main__':unittest.main(verbosity=2)
