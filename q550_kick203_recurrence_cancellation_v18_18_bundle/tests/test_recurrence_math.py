import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from recurrence_math import project,paired_decomposition,instantaneous,train_endpoint_calibrator,ridge_lag

class Tests(unittest.TestCase):
    def setUp(self):
        self.e=np.r_[1.,np.zeros(9)];self.metric=np.ones(10)
        self.ts=np.array([4.5,5.5]);self.lags=np.array([1.,1.5,2.])
    @staticmethod
    def smooth(t):
        t=np.asarray(t,float);p=np.zeros(t.shape+(10,));p[...,1]=np.cos(2.4*t)*(1+.3*np.sin(.9*t));p[...,2]=np.sin(2.4*t)*(.6+.2*np.cos(2.*t));p[...,3]=.2*np.cos(3.1*t);return p
    def test_exact_closure(self):
        d=paired_decomposition(self.smooth,self.ts,self.lags,self.metric,self.e,80)
        self.assertLess(np.max(abs(d['exact']-d['closure'])),1e-12)
    def test_constant_source(self):
        def constant(t):
            t=np.asarray(t);r=np.zeros(t.shape+(10,));r[...,1]=2.;return r
        d=paired_decomposition(constant,self.ts,self.lags,self.metric,self.e,32)
        self.assertTrue(np.allclose(d['exact'],1.,atol=1e-13));self.assertTrue(np.allclose(d['covariance'],0.,atol=1e-13))
    def test_negative_periodic(self):
        def w(t):
            t=np.asarray(t);r=np.zeros(t.shape+(10,));r[...,1]=np.cos(np.pi*t);r[...,2]=np.sin(np.pi*t);return r
        d=paired_decomposition(w,self.ts,np.array([1.]),self.metric,self.e,50)
        self.assertTrue(np.allclose(d['exact'],-1.,atol=1e-13))
    def test_source_weighting_matters(self):
        d=paired_decomposition(self.smooth,self.ts,self.lags,self.metric,self.e,60)
        self.assertGreater(np.max(abs(d['covariance'])),1e-5)
    def test_scramble_changes_age_pairing(self):
        d=paired_decomposition(self.smooth,self.ts,self.lags,self.metric,self.e,64)
        self.assertGreater(np.max(abs(d['scrambled']-d['exact'])),1e-6)
    def test_metric_covariant_with_transform(self):
        m=np.linspace(.6,2.1,10);s=np.sqrt(1/m)
        orig=paired_decomposition(self.smooth,self.ts,self.lags,self.metric,self.e,40)['exact']
        mapped=paired_decomposition(lambda x:self.smooth(x)*s,self.ts,self.lags,m,self.e*s,40)['exact']
        self.assertTrue(np.allclose(orig,mapped,atol=2e-13))
    def test_metric_change_without_transform_does_change_result(self):
        a=paired_decomposition(self.smooth,self.ts,self.lags,self.metric,self.e,40)['exact']
        b=paired_decomposition(self.smooth,self.ts,self.lags,np.linspace(.5,2.,10),self.e,40)['exact']
        self.assertGreater(np.max(abs(a-b)),1e-7)
    def test_projection(self):
        v=self.smooth(self.ts);x=project(v,self.e,self.metric)
        self.assertTrue(np.allclose(x[...,0],0))
    def test_heldout_calibrator_train_only(self):
        e=np.array([[.1,.3],[.2,.4],[.4,.6]]);y=np.array([[.3,.5],[.4,.6],[4.,4.]])
        pred,co=train_endpoint_calibrator(e,y,np.array([4.,6.,8.]),6.)
        self.assertTrue(np.allclose(pred[:2],y[:2]));self.assertGreater(np.max(abs(pred[-1]-y[-1])),1.)
    def test_counterfactual_detects_corrupted_closure(self):
        d=paired_decomposition(self.smooth,self.ts,self.lags,self.metric,self.e,36)
        assert np.max(abs(d['exact']-d['closure']))<1e-12
        broken=d['closure'].copy();broken[0,1]+=.01
        self.assertGreater(np.max(abs(d['exact']-broken)),.009)
    def test_reject_zero_projected_source(self):
        def parallel(t):
            return np.broadcast_to(self.e,np.asarray(t).shape+(10,)).copy()
        with self.assertRaises(ValueError):paired_decomposition(parallel,self.ts,self.lags,self.metric,self.e,20)
    def test_ridge(self):
        a=np.array([[.3,-.6,-.2]])
        self.assertEqual(float(ridge_lag(self.lags,a)[0]),1.5)
if __name__=='__main__':unittest.main()
