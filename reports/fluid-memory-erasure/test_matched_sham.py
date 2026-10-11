import unittest
import numpy as np
from solver import Solver,random_field,curl_gaussian,interventions
from matched_sham import match_translation,norm


class MatchedDisplacementControls(unittest.TestCase):
    def test_equal_change_power_and_translated_pattern(self):
        s=Solver(24,.12);h=random_field(s,10);w=curl_gaussian(s,[np.pi]*3,.45,high_only=True)
        variants,_=interventions(s,np.stack([h+.2*w,h-.2*w]),7,[.13,-.09,.07])
        retained=s.pack(variants['retained']);scr=np.stack([s.pack(variants['scramble'+str(i)]) for i in range(4)])
        moved,matching=match_translation(retained,scr,s.modes_cpu)
        self.assertLess(matching['absolute_matching_error'],1e-12)
        np.testing.assert_allclose(np.sum(abs(moved)**2,axis=1),np.sum(abs(retained)**2,axis=1),atol=1e-17)
        low=np.max(abs(s.modes_cpu),axis=1)<=2
        np.testing.assert_array_equal(moved[:,:,low],retained[:,:,low])
        highret=retained.copy();highret[:,:,low]=0
        highmoved=moved.copy();highmoved[:,:,low]=0
        p=np.random.default_rng(4).uniform(0,2*np.pi,(2,7,3))
        delta=np.array(matching['direction'])*matching['length']
        vm,_=s.sample(s.unpack(highmoved),p,7)
        vr,_=s.sample(s.unpack(highret),p+delta,7)
        np.testing.assert_allclose(vm,vr,atol=1e-13)

if __name__=='__main__':unittest.main(verbosity=2)
