"""Controlled synthetic checks for source recurrence and ridge diagnostics."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ridge_math import project,normalized_pairs,fixed_support,instantaneous,variable_support,ridge_row

def source_v(arr,dim=10,omega=0):
    x=np.asarray(arr,float)
    ans=np.zeros(x.shape+(dim,))
    ans[...,1]=np.cos(omega*x) if omega else 1.
    return ans

class TestRidgeMath(unittest.TestCase):
    def setUp(self):self.w=np.ones(10);self.e=np.r_[1.,np.zeros(9)]
    def test_constant_source_matches_any_lag(self):
        k=fixed_support(source_v,np.array([4.,6.]),np.array([1.,1.5,2.]),60,self.w,self.e)
        self.assertTrue(np.allclose(k,1.,atol=1e-13))
    def test_sinusoidal_antagonism_at_one_xi(self):
        g=lambda x:source_v(x,omega=np.pi)
        k=fixed_support(g,np.array([4.,6.]),np.array([1.]),60,self.w,self.e)
        self.assertTrue(np.allclose(k,-1.,atol=1e-12))
    def test_metric_covariance_with_consistent_coordinates(self):
        def g(x):
            x=np.asarray(x);r=source_v(x);r[...,2]=np.sin(x);r[...,3]=.7*np.cos(1.7*x);return r
        x=np.array([4.0,5.5]);lags=np.array([1.0,1.5,2.0]);b=np.linspace(.7,2.,10)
        k1=fixed_support(g,x,lags,80,np.ones(10),self.e)
        k2=fixed_support(lambda t:g(t)*b,x,lags,80,b**-2,self.e*b)
        self.assertTrue(np.allclose(k1,k2,atol=1e-12))
    def test_project_is_transverse(self):
        w=np.linspace(.5,2.,10);p=project(self.e,np.ones((4,10)),w)
        self.assertLess(np.max(abs(p@(w*self.e))),1e-12)
    def test_instantaneous_source_returns_negative_one(self):
        g=lambda x:source_v(x,omega=np.pi)
        k=instantaneous(g,np.array([4.,5.]),np.array([1.]),self.w,self.e)
        self.assertTrue(np.allclose(k,-1.))
    def test_variable_support_rejects_insufficient_history(self):
        k=variable_support(source_v,np.array([3.,4.]),np.array([2.]),self.w,self.e)
        self.assertTrue(np.isnan(k[0,0]));self.assertAlmostEqual(k[1,0],1.,places=12)
    def test_ridge_parabolic_interpolation(self):
        lags=np.linspace(1.,2.,41)
        z=ridge_row(lags,(lags-1.383)**2-.9)
        self.assertAlmostEqual(z['lag_quad'],1.383,places=10)
        self.assertFalse(z['edge'])
    def test_ridge_flags_boundary(self):
        x=np.linspace(1.,2.,21);r=ridge_row(x,x)
        self.assertTrue(r['edge']);self.assertEqual(r['lag_quad'],1.)
    def test_midpoint_convergence_smooth(self):
        def g(x):
            x=np.asarray(x);r=source_v(x,omega=2.2);r[...,2]=.2*np.sin(3*x);return r
        xs=np.array([4.25,6.75]);Ls=np.array([1.2,1.5])
        k1=fixed_support(g,xs,Ls,80,self.w,self.e)
        k2=fixed_support(g,xs,Ls,160,self.w,self.e)
        self.assertLess(np.max(abs(k1-k2)),.001)

if __name__=='__main__':unittest.main()
