import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from chirality import signed_features


def sample():
    rng=np.random.default_rng(602)
    start=rng.uniform(-.4,.4,(32,3))
    return np.array([start+t*rng.normal(0,.05,start.shape) for t in np.linspace(0,1,9)])


def test_signed_reflection_rotation_and_reversal():
    x=sample(); c,h=signed_features(x)
    assert abs(c).max()>1e-5 and abs(h[:,0]).max()>1e-5
    R=np.diag([-1,1,1])
    cm,hm=signed_features(x@R)
    np.testing.assert_allclose(cm,-c,atol=1e-12)
    np.testing.assert_allclose(hm,h*[-1,1,1],atol=1e-12)
    cr,hr=signed_features(x,reverse=True)
    np.testing.assert_allclose(cr,c,atol=1e-12)
    np.testing.assert_allclose(hr,h*[-1,1,-1],atol=1e-12)
    Q=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    cq,hq=signed_features(x@Q)
    np.testing.assert_allclose(cq,c,atol=1e-12)
    np.testing.assert_allclose(hq,h,atol=1e-12)


def test_identity_translation_and_causality():
    x=sample(); c,h=signed_features(x)
    perm=np.random.default_rng(603).permutation(x.shape[1])
    cp,hp=signed_features(x[:,perm]+[.2,-.5,6.1])
    np.testing.assert_allclose(cp,c[perm],atol=1e-11)
    np.testing.assert_allclose(hp,h[perm],atol=1e-11)
    full=np.concatenate([x,np.full_like(x,np.nan)])
    cf,hf=signed_features(full[:len(x)])
    np.testing.assert_array_equal(cf,c)
    np.testing.assert_array_equal(hf,h)


def test_planar_and_static_histories():
    x=sample(); x[:,:,2]=0
    c,h=signed_features(x)
    np.testing.assert_allclose(c,0,atol=1e-13)
    np.testing.assert_allclose(h,0,atol=1e-13)
    _,h=signed_features(np.repeat(sample()[:1],9,axis=0))
    np.testing.assert_allclose(h,0,atol=1e-13)
    with pytest.raises(ValueError): signed_features(sample()[:8])
