import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from raw_jet_math import *

class RawSampleJetTest(unittest.TestCase):
    def test_basis_endpoints(self):
        t=basis_times(.2,33)
        self.assertAlmostEqual(t[0],-.2)
        self.assertAlmostEqual(t[-1],.2)
        self.assertEqual(len(t),33)
    def test_invalid_even_count(self):
        with self.assertRaises(ValueError):basis_times(.1,14)
    def test_quintic_polynomial_exact_recovery(self):
        t=basis_times(.5,33);_,B=polynomial_operator(t,5,.5)
        y=sum((t/.5)**k*(k+1) for k in range(6))
        np.testing.assert_allclose(B@y,np.arange(1,7),rtol=1e-12,atol=1e-12)
    def test_polynomial_coefficient_covariance(self):
        t=basis_times(.5,33);_,B=polynomial_operator(t,5,.5)
        K=B[:4]@B[:4].T
        self.assertTrue(np.min(np.linalg.eigvalsh(K))>0)
        self.assertGreater(abs(K[0,2]),.001)
    def test_whitened_covariance(self):
        t=basis_times(.2,33)
        jets=[np.ones((3,10)) for _ in range(4)]
        b,c,w,o,h=jet_compression(t,.2,5,3,jets)
        self.assertLess(test_covariance(c,3,w),1e-10)
    def test_compressed_polynomial_signal(self):
        t=basis_times(.2,33)
        H=np.eye(10)[:3]
        jets=[H,np.zeros_like(H),np.zeros_like(H),np.zeros_like(H)]
        b,c,w,o,h=jet_compression(t,.2,5,3,jets)
        v=np.arange(10)/10
        z=np.tile(H@v,(33,1)).reshape(-1)
        np.testing.assert_allclose(compressed_observation(z,3,b,w),o@v,rtol=1e-12,atol=1e-12)
    def test_direct_constant_map(self):
        P=np.stack([np.eye(10)]*5)
        H=np.eye(10)[:4]
        O=direct_observation(P,2,[0,1,2,3,4],H)
        np.testing.assert_allclose(O,np.tile(H,(5,1)))
    def test_raw_least_squares_exact(self):
        rng=np.random.default_rng(3);O=rng.normal(size=(33,10));x=rng.normal(size=10)
        np.testing.assert_allclose(direct_estimator(O)@(O@x),x,rtol=1e-12,atol=1e-12)
    def test_correlated_coefficient_noise_monte_carlo(self):
        t=basis_times(.2,33);_,B=polynomial_operator(t,5,.2)
        cov=B[:4]@B[:4].T
        x=np.random.default_rng(32).normal(size=(4000,33))
        empirical=np.cov(x@B[:4].T,rowvar=False)
        self.assertLess(np.linalg.norm(empirical-cov)/np.linalg.norm(cov),.08)
    def test_nontrivial_model_truncation(self):
        t=basis_times(.3,33)
        V,B=polynomial_operator(t,3,.3)
        y=t**4
        self.assertGreater(np.linalg.norm(B[:4]@y),1e-6)
    def test_continuous_time_no_node_duplication(self):
        from scipy.interpolate import CubicSpline
        grid=np.linspace(3.8,4.2,33)
        P=np.stack([(1.+.2*(t-4.))*np.eye(10) for t in grid])
        spline=CubicSpline(grid,P,axis=0)
        H=np.eye(10)[:3]
        times=4.+basis_times(.05,33)
        O=continuous_observation(spline,np.eye(10),times,H)
        self.assertEqual(len(np.unique(times)),33)
        self.assertEqual(len(np.unique(O[::3,0])),33)
        np.testing.assert_allclose(O[::3,:],
          np.array([(1.+.2*(t-4.))*H[0] for t in times]),atol=1e-12)
    def test_unknown_jet_covariance_not_assumed_identity(self):
        t=basis_times(.2,33);_,B=polynomial_operator(t,5,.2)
        C=B[:4]@B[:4].T
        self.assertGreater(np.linalg.norm(C-np.eye(4)),.1)

if __name__=='__main__':unittest.main()
