import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from noise_gauge_math import oracle_stats,gauge_test,orthogonal,gauge_matrices

class TestNoiseGauge(unittest.TestCase):
    def test_unitary(self):
        Q=orthogonal(10);self.assertLess(np.linalg.norm(Q.T@Q-np.eye(10)),1e-13)
    def test_gauge_prediction(self):
        rng=np.random.default_rng(1); O=rng.normal(size=(15,10))
        R=gauge_test(O,3,4)
        self.assertLess(R['physical_prediction_error'],1e-10)
        self.assertGreater(R['cond_T'],1.)
    def test_source_gauge_not_automatically_orthogonal(self):
        T,B=gauge_matrices(10,3,3)
        self.assertGreater(np.linalg.cond(T),50.)
    def test_observer_transform_invertible(self):
        T,B=gauge_matrices(10,3,3)
        self.assertEqual(np.linalg.matrix_rank(T),10)
        self.assertEqual(np.linalg.matrix_rank(B),12)
    def test_noise_analytic(self):
        O=np.diag(np.linspace(1.,2.,10))
        out=oracle_stats(O,[.01],seed=1,trials=10000)[0]
        self.assertLess(abs(out['empirical_rmse']/out['analytical_expected_rmse']-1),.02)
    def test_noise_snr_monotone(self):
        O=np.diag(np.linspace(.1,2.,10))
        out=oracle_stats(O,[1e-6,1e-3],seed=1)
        self.assertGreater(out[0]['weak_direction_snr'],out[1]['weak_direction_snr'])
    def test_covariance_consistent_rescaling(self):
        rng=np.random.default_rng(2);O=rng.normal(size=(12,10));m=3
        sig=np.repeat([1e-4,2e-4,5e-3,1e-2],m)
        s0=np.linalg.svd(O/sig[:,None],compute_uv=False)
        for tau in [.03,.1,.2]:
            d=np.repeat(tau**np.arange(4),m)
            s=np.linalg.svd((O*d[:,None])/(sig*d)[:,None],compute_uv=False)
            self.assertLess(np.max(abs(s-s0)),1e-9)
    def test_whitened_gauge(self):
        rng=np.random.default_rng(3);O=rng.normal(size=(12,10))
        T,B=gauge_matrices(10,3,3)
        s1=np.linalg.svd(O,compute_uv=False)
        O2=np.linalg.solve(B,B@O@np.linalg.inv(T))@T
        s2=np.linalg.svd(O2,compute_uv=False)
        self.assertLess(np.linalg.norm(s1-s2),1e-9)
    def test_jet_scale_not_new_information(self):
        a=np.array([1.,2.,3.]);s=np.array([.1,.2,.5]);d=np.array([.01,.1,1.])
        self.assertTrue(np.allclose(d*a/(d*s),a/s))
    def test_covariance_positive(self):
        T,B=gauge_matrices(10,3,3);eig=np.linalg.eigvalsh(B@B.T)
        self.assertGreater(eig.min(),0.)

if __name__=='__main__':unittest.main()
