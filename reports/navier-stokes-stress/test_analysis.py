"""Analytic tests for full-spectrum comparisons and first-failure censoring."""
import json, unittest
from pathlib import Path
import numpy as np
from solver import Solver
from analyze import relative,prefix

class AnalysisTests(unittest.TestCase):
    def test_common_and_missing_modes_count_once(self):
        coarse=Solver(24,.01);fine=Solver(48,.01)
        a=coarse.initial('abc');base=fine.initial('abc')
        x=np.arange(48)*fine.dx
        u=np.zeros((3,48,48,48));u[1]=np.sin(10*x)[:,None,None]
        high=fine.project(fine.hat(u))
        self.assertLess(relative(a,base,coarse),1e-14)
        for scale in [1.,1.1]:
            b=scale*base+.2*high
            expected=np.sqrt(((scale-1)**2*fine.inner(base,base)+.04*fine.inner(high,high))/fine.inner(b,b))
            self.assertAlmostEqual(relative(a,b,coarse),expected,places=13)
        self.assertAlmostEqual(relative(a*.9,a,coarse),.1,places=13)

    def test_later_pass_cannot_erase_failure(self):
        self.assertEqual(prefix([dict(t=0,pass_=True),dict(t=.5,pass_=False),dict(t=1,pass_=True)],'pass_'),0)
        self.assertIsNone(prefix([dict(t=0,pass_=False)],'pass_'))

if __name__=='__main__':
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AnalysisTests))
    Path('analysis-validation.json').write_text(json.dumps(dict(tests=r.testsRun,passed=r.wasSuccessful(),failures=len(r.failures),errors=len(r.errors)),indent=2))
    raise SystemExit(not r.wasSuccessful())
