import unittest,numpy as np
from corrected_solver import Solver,interventions,random_field
from final_experiment import basis
class IntegerWaveControls(unittest.TestCase):
    def test_grid_independent_reset_at_degenerate_modes(self):
        values=[]
        # The previous n=80 implementation failed at k=(2,4,6), where
        # the primary cross product vanishes analytically.
        for n in [64,80,96]:
            s=Solver(n,.02);h=random_field(s,80000);q=basis(s,[np.pi]*3,.45,4,True)[2]
            v,_=interventions(s,np.stack([h+.2*q,h-.2*q]),80000,[.13,-.09,.07],4)
            modes=s.modes_cpu;coeff=s.pack(v['reset'])
            select=[tuple(k) for k in modes].index((2,4,6))
            values.append(coeff[:,:,select])
            np.testing.assert_array_equal(s.k,s.k_cpu)
            self.assertLess(abs(np.sum(s.k*v['reset'],axis=1)).max()/n**3,1e-14)
        for a in values[1:]:np.testing.assert_allclose(a,values[0],atol=2e-15,rtol=1e-12)
    def test_reset_same_polynomial_with_explicit_modal_power(self):
        values=[]
        for n in [64,80,96]:
            s=Solver(n,.02);h=np.zeros((2,)+s.shape,complex)
            for k in [(2,4,6),(-2,-4,6),(3,6,9)]:
                j=np.where(np.all(s.modes_cpu==k,axis=1))[0][0]
                c=np.zeros((2,3,len(s.modes_cpu)),complex);c[:,0,j]=.01;c[:,1,j]=-.005
                h+=s.project(s.unpack(c))
            v,_=interventions(s,h,1,[.13,-.09,.07],4)
            coeff=s.pack(v['reset']);lookup={tuple(k):j for j,k in enumerate(s.modes_cpu)}
            values.append(np.stack([coeff[:,:,lookup[k]] for k in [(2,4,6),(-2,-4,6),(3,6,9)]]))
        for a in values[1:]:np.testing.assert_allclose(a,values[0],atol=2e-15)
if __name__=='__main__':unittest.main(verbosity=2)
