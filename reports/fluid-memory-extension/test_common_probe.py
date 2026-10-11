import unittest,numpy as np
from experiment import DrivenSolver,basis
from reference_solver import random_field,interventions
from common_probe import construct
class CommonProbeControl(unittest.TestCase):
    def test_same_pulse_same_strength_same_work_all_arms(self):
        s=DrivenSolver(24,.12,backend='cpu');h=random_field(s,81000);writer=basis(s,[np.pi]*3,.45,4,True)[2]
        v,_=interventions(s,np.stack([h+.2*writer,h-.2*writer]),81000,[.13,-.09,.07],4)
        states=np.concatenate(list(v.values()));pulse,audit=construct(s,states,[np.pi]*3)
        self.assertLess(audit['work_error'],1e-12);self.assertLess(audit['rms_error'],1e-14)
        actual=.5*s.inner(states+pulse,states+pulse)-.5*s.inner(states,states)
        np.testing.assert_allclose(actual,.0002,atol=1e-13)
        self.assertLess(abs(np.sum(s.k*pulse,axis=0)).max()/s.n**3,1e-13)
if __name__=='__main__':unittest.main(verbosity=2)
