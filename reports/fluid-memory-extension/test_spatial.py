import unittest,numpy as np
from spatial_assessment import footprint
class SpatialControls(unittest.TestCase):
    def test_known_shear_rms_and_radial_cumulative_bounds(self):
        modes=np.array([[0,1,0],[0,-1,0]])
        coeff=np.zeros((1,3,2),complex);coeff[0,0]=[-.5j,.5j]
        r=footprint(coeff,modes,[np.pi]*3)[0]
        self.assertAlmostEqual(r['rms'],np.sqrt(.5),places=13)
        self.assertTrue(all(0<=x<=1 for x in r['radial_fraction']))
        self.assertTrue(np.all(np.diff(r['radial_fraction'])>=0))
        self.assertAlmostEqual(r['outside_radius2'],1-r['radial_fraction'][3])
    def test_zero_response_has_no_defined_spatial_extent(self):
        r=footprint(np.zeros((1,3,1),complex),np.array([[0,0,0]]),[np.pi]*3)[0]
        self.assertEqual(r['rms'],0);self.assertIsNone(r['outside_radius2'])
if __name__=='__main__':unittest.main(verbosity=2)
