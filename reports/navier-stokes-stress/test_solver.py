import json, unittest
from pathlib import Path
import numpy as np
from solver import Solver,restrict,difference

class ScienceTests(unittest.TestCase):
    def test_incompressible_and_energy_budget(self):
        for family in ['kida','taylor_green','abc']:
            s=Solver(24,.005);h=s.initial(family,.001,7);r,d,_=s.rhs(h)
            self.assertLess(abs(s.inner(h,r)+d),1e-12)
            div=1j*sum(k*c for k,c in zip(s.k,h));self.assertLess(s.inner(div,div),1e-25)
            self.assertLess(np.max(abs(h[:,~s.keep])),1e-12)
    def test_same_start_and_energy_across_grids(self):
        for family in ['kida','taylor_green']:
            a=Solver(24,.01);b=Solver(32,.01);h=a.initial(family,.01,3);g=b.initial(family,.01,3)
            self.assertLess(np.max(abs(restrict(g,24)-h))/24**3,1e-14)
            self.assertAlmostEqual(a.inner(h,h),a.inner(a.initial(family),a.initial(family)),places=13)
            self.assertLess(difference(h,g,a),3e-8)
    def test_exact_decay(self):
        for family,k2 in [('abc',1),('taylor_green_2d',2)]:
            s=Solver(16,.1);h=s.initial(family);start=h.copy();loss=0
            for _ in range(20):h,dl,_=s.step(h,.025);loss+=dl
            target=start*np.exp(-.1*k2*.5)
            self.assertLess(np.sqrt(s.inner(h-target,h-target)/s.inner(target,target)),1e-11)
            self.assertLess(abs(.5*s.inner(h,h)-.5*s.inner(start,start)+loss),1e-11)
    def test_fourth_order_time_and_translation(self):
        s=Solver(16,.05);start=s.initial('abc');mean=np.array([1.,.7,.3]);start[:,0,0,0]=mean*16**3
        target=(start*np.exp(-.05*s.k2-.0j)*np.exp(-1j*sum(k*v for k,v in zip(s.k,mean))))
        target[:,0,0,0]=start[:,0,0,0];errors=[]
        for dt in [.2,.1,.05]:
            h=start.copy()
            for _ in range(round(1/dt)):h,_,_=s.step(h,dt)
            errors.append(np.sqrt(s.inner(h-target,h-target)))
        self.assertGreater(errors[0]/errors[1],14);self.assertGreater(errors[1]/errors[2],14)
    def test_nonpositive_viscosity_rejected(self):
        with self.assertRaises(ValueError):Solver(16,0)

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ScienceTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    Path('validation.json').write_text(json.dumps(dict(tests=result.testsRun,passed=result.wasSuccessful(),failures=len(result.failures),errors=len(result.errors)),indent=2))
    raise SystemExit(not result.wasSuccessful())
