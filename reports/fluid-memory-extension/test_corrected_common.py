import unittest,numpy as np
from final_experiment import DrivenSolver,basis
from corrected_solver import random_field,interventions
from common_probe import construct
class CommonProbeControl(unittest.TestCase):
    def test_state_label_permutation_preserves_the_common_impulse(self):
        s=DrivenSolver(24,.12,backend='cpu');h=random_field(s,81000);writer=basis(s,[np.pi]*3,.45,4,True)[2]
        v,_=interventions(s,np.stack([h+.2*writer,h-.2*writer]),81000,[.13,-.09,.07],4)
        states=np.concatenate(list(v.values()));pulse,a=construct(s,states,[np.pi]*3)
        permutation=np.random.default_rng(29).permutation(len(states));other,b=construct(s,states[permutation],[np.pi]*3)
        np.testing.assert_allclose(pulse/s.n**3,other/s.n**3,atol=1e-12)
        self.assertEqual(a['constraint_rank'],b['constraint_rank'])
    def test_degenerate_preference_fallback_is_label_independent(self):
        s=DrivenSolver(24,.12,backend='cpu');q=basis(s,[np.pi]*3)[2]
        states=np.stack([.1*q,.2*q,-.3*q]);pulse,a=construct(s,states,[np.pi]*3)
        other,b=construct(s,states[[2,0,1]],[np.pi]*3)
        self.assertTrue(a['fallback']);self.assertTrue(b['fallback'])
        self.assertEqual(a['fallback_basis_index'],b['fallback_basis_index'])
        np.testing.assert_allclose(pulse/s.n**3,other/s.n**3,atol=1e-12)
        self.assertLess(a['work_error'],1e-12)
    def test_same_pulse_same_strength_same_work_all_arms(self):
        s=DrivenSolver(24,.12,backend='cpu');h=random_field(s,81000);writer=basis(s,[np.pi]*3,.45,4,True)[2]
        v,_=interventions(s,np.stack([h+.2*writer,h-.2*writer]),81000,[.13,-.09,.07],4)
        states=np.concatenate(list(v.values()));pulse,audit=construct(s,states,[np.pi]*3)
        self.assertLess(audit['work_error'],1e-12);self.assertLess(audit['rms_error'],1e-14)
        actual=.5*s.inner(states+pulse,states+pulse)-.5*s.inner(states,states)
        np.testing.assert_allclose(actual,.0002,atol=1e-13)
        self.assertLess(abs(np.sum(s.k*pulse,axis=0)).max()/s.n**3,1e-13)
if __name__=='__main__':unittest.main(verbosity=2)
