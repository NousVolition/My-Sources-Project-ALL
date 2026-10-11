import unittest
import numpy as np
from analyze import spectral_distance,features,interval,P
from run import marker_cloud


class AnalysisControls(unittest.TestCase):
    def test_missing_fourier_mode_counted(self):
        k=np.array([[0,0,0],[1,0,0]]);kf=np.array([[0,0,0],[1,0,0],[0,0,1]])
        a=np.zeros((1,3,2),complex);b=np.zeros((1,3,3),complex)
        a[0,0,1]=b[0,0,1]=1;b[0,0,2]=1
        self.assertAlmostEqual(spectral_distance(a,k,b,kf)[0],np.sqrt(2/3))

    def test_field_error_near_roundoff(self):
        k=np.array([[1,0,0]]);a=np.ones((1,3,1),complex);b=a*(1+1e-12)
        self.assertAlmostEqual(spectral_distance(a,k,b,k)[0],1e-12,delta=2e-16)

    def test_static_history_is_zero_and_identity_independent(self):
        x=marker_cloud(1);times=np.arange(5)*.1;pos=np.tile(x,(5,1,1));v=np.zeros_like(x);g=np.zeros((10,3,3))
        a,h=features(times,pos,v,g)
        # Cosine summaries equal one for static history; all rate features zero.
        self.assertLess(abs(h.reshape(10,3,-1)[...,:4]).max(),1e-14)
        order=np.array([7,2,9,0,1,3,4,5,6,8])
        probes=np.concatenate([np.arange(10+6*i,10+6*i+6) for i in order])
        reindex=np.concatenate([order,probes])
        ar,hr=features(times,pos[:,reindex],v[reindex],g[order])
        np.testing.assert_allclose(ar,a[order],atol=1e-14);np.testing.assert_allclose(hr,h[order],atol=1e-14)

    def test_future_access_rejected(self):
        x=marker_cloud(1)
        with self.assertRaises(ValueError):features(np.arange(6)*.1,np.tile(x,(6,1,1)),np.zeros_like(x),np.zeros((10,3,3)))
        with self.assertRaises(ValueError):features(np.array([0,.1,.5,.2,.4]),np.tile(x,(5,1,1)),np.zeros_like(x),np.zeros((10,3,3)))

    def test_split_seeds_disjoint(self):
        all_sets=[set(P[k]) for k in ['train_seeds','validation_seeds','test_seeds','pilot_seeds']]
        for i,a in enumerate(all_sets):
            for b in all_sets[i+1:]:self.assertFalse(a&b)

    def test_bootstrap_uses_flows(self):
        r=interval([1,2,3],reps=100)
        self.assertEqual(r['n_independent_flows'],3);self.assertEqual(r['mean'],2)

if __name__=='__main__':unittest.main(verbosity=2)
