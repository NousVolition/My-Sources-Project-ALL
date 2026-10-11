import unittest,numpy as np
from followup_analysis import claim,common_measure,common_response
from analysis import distance
from final_experiment import P
class InferenceControls(unittest.TestCase):
    def test_response_invariant_to_marker_identity_and_history_names(self):
        # Known diagonal deformation maps; paired probe minus base rates.
        base=np.tile(np.eye(3),(2,14,10,1,1));kicked=base.copy()
        for j in range(14):
            rates=np.linspace(.001,.004,10)*(1+j%2)
            kicked[1,j,:,0,0]=np.exp(.4*rates);kicked[1,j,:,1,1]=np.exp(-.4*rates)
        a=dict(time=np.array([0.,.4]),tangent=kicked,base_tangent=base)
        before=common_measure(a,{})
        perm=np.random.default_rng(4).permutation(10)
        b=dict(a,tangent=kicked[:,:,perm],base_tangent=base[:,:,perm])
        np.testing.assert_allclose(list(before.values()),list(common_measure(b,{}).values()),atol=1e-15)
        hist=np.arange(14).reshape(7,2)[:,::-1].reshape(-1)
        b=dict(a,tangent=kicked[:,hist],base_tangent=base[:,hist])
        np.testing.assert_allclose(list(before.values()),list(common_measure(b,{}).values()),atol=1e-15)
        np.testing.assert_allclose(common_response(a,{})['retained'][1]-common_response(a,{})['retained'][0],np.linspace(.001,.004,10),atol=1e-14)
    def test_identical_full_states_have_zero_history_contrast(self):
        f=np.tile(np.eye(3),(2,14,10,1,1));g=f.copy();g[1,:,:,0,0]=np.exp(.4*.002)
        a=dict(time=np.array([0.,.4]),tangent=g,base_tangent=f)
        self.assertTrue(all(x==0 for x in common_measure(a,{}).values()))
    def test_statistical_signal_requires_resolution_before_physical_claim(self):
        low=dict(mean=1.,ci95=[.8,1.2]);high=dict(mean=1.,ci95=[.8,1.2])
        self.assertFalse(claim(low,high,[.2,.2])['effect_resolution_passed'])
        self.assertTrue(claim(low,high,[.01,-.01])['effect_resolution_passed'])
        high=dict(mean=1.,ci95=[-.1,1.2]);self.assertFalse(claim(low,high,[0.,0.])['effect_resolution_passed'])
    def test_spectral_distance_accounts_for_unresolved_modes(self):
        ka=np.array([[0,0,0]]);kb=np.array([[0,0,0],[0,0,1]])
        a=np.array([[1.+0j]]);b=np.array([[1.+0j,1.+0j]])
        self.assertAlmostEqual(float(distance(a,ka,b,kb)),np.sqrt(2/3))
if __name__=='__main__':unittest.main(verbosity=2)
