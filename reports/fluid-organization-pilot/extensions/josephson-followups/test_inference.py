import numpy as np
from prediction import fit_kernel,predict,shuffled_history
from ramps import trajectory

def test_kernel_recovers_smooth_holdout():
    x=np.linspace(-1,1,30)[:,None];v=np.linspace(-.95,.95,12)[:,None]
    m=fit_kernel(x,v,x[:,0]**2,v[:,0]**2)
    assert np.sqrt(np.mean((predict(m,v)-v[:,0]**2)**2))<.01
def test_kernel_transform_only_uses_training():
    x=np.arange(10.)[:,None];v=np.array([[1000.],[2000.]])
    m=fit_kernel(x,v,np.sin(x[:,0]),np.zeros(2))
    np.testing.assert_allclose(m['mean'],x.mean(0));np.testing.assert_allclose(m['scale'],x.std(0))
def test_shuffle_preserves_run_slots_and_strata():
    d={'groups':np.repeat(np.arange(8),3),'history':np.arange(24).reshape(-1,1),'stratum':np.repeat(np.arange(8)//4,3)}
    q=shuffled_history(d,np.random.default_rng(1)).reshape(8,3)
    np.testing.assert_array_equal(np.diff(q,axis=1),np.ones((8,2)))
    np.testing.assert_array_equal(q[:,0]//12,np.arange(8)//4)
    assert np.all(q[:,0]!=np.arange(8)*3)
def test_staircase_handles_steps_smaller_than_sample_bins():
    s=np.array([[.1,-1.]])
    t,b,q=trajectory(s,2,4,.9,1.3,.5,step=.0125)
    assert q.shape==(9,1,2) and np.isfinite(q).all()
