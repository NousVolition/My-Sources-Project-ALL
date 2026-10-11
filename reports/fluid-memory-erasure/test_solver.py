"""Independent analytic and implementation controls; run before production."""
import unittest
import numpy as np
from solver import Solver, random_field, curl_gaussian, phase_values, interventions


class ScientificControls(unittest.TestCase):
    def test_exact_spectral_velocity_and_gradient(self):
        s=Solver(24,.12)
        a=np.arange(s.n)*s.dx
        x,y,z=np.meshgrid(a,a,a,indexing='ij')
        u=np.array([np.sin(y),np.cos(z),np.sin(x)])
        h=s.project(s.hat(u)[None])
        p=np.random.default_rng(9).uniform(0,2*np.pi,(1,13,3))
        v,g=s.sample(h,p,13)
        true=np.stack([np.sin(p[...,1]),np.cos(p[...,2]),np.sin(p[...,0])],axis=-1)
        exact=np.zeros_like(g);exact[...,0,1]=np.cos(p[...,1]);exact[...,1,2]=-np.sin(p[...,2]);exact[...,2,0]=np.cos(p[...,0])
        np.testing.assert_allclose(v,true,atol=3e-14)
        np.testing.assert_allclose(g,exact,atol=3e-14)
        self.assertLess(abs(np.trace(g,axis1=-2,axis2=-1)).max(),1e-14)

    def shear_run(self,dt):
        s=Solver(16,.2)
        axes=np.arange(s.n)*s.dx
        x,y,z=np.meshgrid(axes,axes,axes,indexing='ij')
        u=np.zeros((3,s.n,s.n,s.n));u[0]=np.sin(y)
        h=s.hat(u)[None];p=np.array([[[.3,.7,.5],[.8,1.2,1.]]]);initial=p.copy()
        f=np.tile(np.eye(3),(1,2,1,1));budget=np.zeros(1)
        for _ in range(round(.4/dt)):
            h,p,f,db,_=s.step(h,p,f,dt);budget+=db
        factor=-np.expm1(-.2*.4)/.2
        true=initial.copy();true[...,0]+=np.sin(initial[...,1])*factor
        exact=np.tile(np.eye(3),(1,2,1,1));exact[...,0,1]=np.cos(initial[...,1])*factor
        return max(abs(p-true).max(),abs(f-exact).max()), h,s,budget

    def test_coupled_shear_solution_and_fourth_order(self):
        errors=[self.shear_run(dt)[0] for dt in [.1,.05,.025]]
        self.assertGreater(errors[0]/errors[1],14)
        self.assertGreater(errors[1]/errors[2],14)
        self.assertLess(errors[-1],1e-11)

    def test_energy_parseval_and_budget(self):
        _,h,s,b=self.shear_run(.025)
        e=.5*s.inner(h,h)[0]
        self.assertAlmostEqual(e,float(.5*np.mean(np.sum(s.real(h)[0]**2,axis=0))),places=14)
        self.assertLess(abs(e-.25-b[0]),1e-10)

    def test_grid_independent_initial_field(self):
        p=np.random.default_rng(3).uniform(0,2*np.pi,(1,9,3))
        observed=[]
        for n in [24,32,48]:
            s=Solver(n,.12);h=random_field(s,60000)[None];observed.append(s.sample(h,p,9))
        for v,g in observed[1:]:
            np.testing.assert_allclose(v,observed[0][0],atol=1e-14)
            np.testing.assert_allclose(g,observed[0][1],atol=5e-14)

    def test_phase_reality_and_grid_independence(self):
        k=np.array([[1,2,0],[-1,-2,0],[0,-3,1],[0,0,0]])
        phases=phase_values(k,1)
        self.assertAlmostEqual(phases[0],-phases[1])
        np.testing.assert_array_equal(phase_values(k[[2,0,1]],1),phases[[2,0,1]])
        s=Solver(24,.12);h=random_field(s,4)
        d=curl_gaussian(s,[np.pi]*3,.45,high_only=True)
        states,checks=interventions(s,np.stack([h+.2*d,h-.2*d]),71,[.13,-.09,.07])
        for name,a in states.items():
            np.testing.assert_allclose(s.hat(s.real(a)),a,atol=1e-11)
            self.assertLess(checks[name]['power_error'],1e-14)
            self.assertLess(checks[name]['coarse_error'],1e-14)
            self.assertLess(checks[name]['divergence_error'],1e-13)
        np.testing.assert_array_equal(states['reset'][0],states['reset'][1])

    def test_pack_roundtrip(self):
        s=Solver(24,.12);h=random_field(s,5)[None]
        np.testing.assert_allclose(s.unpack(s.pack(h)),h,atol=1e-12)

    def test_no_probe_identical_branch(self):
        s=Solver(16,.12);single=random_field(s,6);h=np.stack([single,single])
        x=np.tile(np.array([[[.7,.4,.8]]]),(2,1,1));f=np.tile(np.eye(3),(2,1,1,1))
        h,x,f,_,_=s.step(h,x,f,.01)
        np.testing.assert_array_equal(h[0],h[1]);np.testing.assert_array_equal(x[0],x[1]);np.testing.assert_array_equal(f[0],f[1])

    def test_positive_viscosity(self):
        with self.assertRaises(ValueError):Solver(24,0)


if __name__=='__main__':
    unittest.main(verbosity=2)
