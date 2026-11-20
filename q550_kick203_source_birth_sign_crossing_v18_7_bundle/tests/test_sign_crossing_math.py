"""Synthetic tests only; do not interpret as Q550 results."""
import unittest,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from sign_crossing_math import (cumulative_simpson,continuous_integral,
 sign_and_metric_rates,roots_from_grid,polar_rotation,kinematic_rate_check,angle_deg)

class TestSignCrossing(unittest.TestCase):
    def test_simpson_polynomial(self):
        x=np.linspace(0,8,641);g=np.zeros((len(x),11));g[:,0]=3*x*x+2*x+1
        I=cumulative_simpson(x,g,2.5)
        self.assertTrue(np.allclose(I[:,0],2.5*(x**3+x*x+x),rtol=0,atol=3e-12))
    def test_hermite_cubic(self):
        x=np.linspace(0,8,641);g=np.zeros((len(x),11));g[:,0]=3*x*x+2*x+1
        I=cumulative_simpson(x,g,1.)
        p=continuous_integral(x,g,I,1.)
        for y in [0.49,2.371,7.919]:self.assertAlmostEqual(p(y)[0],y**3+y*y+y,places=10)
    def test_roots_all(self):
        x=np.linspace(0,8,641)
        roots,_=roots_from_grid(x,lambda z:(z-3.37)*(z-6.61),start=1)
        self.assertEqual(len(roots),2)
        self.assertTrue(np.allclose(roots,[3.37,6.61],atol=1e-10))
    def test_rates_nontrivial_strain(self):
        E=np.eye(11)[:,0]; T=np.zeros(11);T[:2]=[1,0.6]
        J=np.zeros((11,11));J[0,0]=0.3;J[1,1]=-0.4
        P=np.eye(11);P[0,0]=1.2;P[1,1]=0.8
        g=np.zeros(11);g[:2]=[-.5,.25]
        _,metric,birth=sign_and_metric_rates(P,J,2.5,E,T,g)
        dP=2.5*J@P
        dq=2.5*(J@(P@T)+P@g)
        de=2.5*J@(P@E)
        expected=de@(P@T)+(P@E)@dq
        self.assertAlmostEqual(metric+birth,expected,places=12)
        self.assertGreater(abs(metric),0.1)
        self.assertGreater(abs(birth),0.1)
    def test_rotation_and_11d_gauge(self):
        A=np.random.default_rng(7).normal(size=(11,11)); U,_,Vt=np.linalg.svd(A)
        R=U@Vt
        P=R@np.diag(np.linspace(.5,1.5,11))
        self.assertLess(np.linalg.norm(polar_rotation(P)-R),1e-12)
        self.assertLess(np.linalg.norm(polar_rotation(P).T@polar_rotation(P)-np.eye(11)),1e-12)
        E=np.arange(1.,12.);T=np.arange(11.,0.,-1)
        self.assertAlmostEqual((P@E)@(P@T),(R.T@P@E)@(R.T@P@T),places=10)
    def test_co_moving_returns_constant(self):
        P=np.diag(np.arange(1.,12.));E=np.arange(1.,12.)
        self.assertLess(np.linalg.norm(np.linalg.solve(P,P@E)-E),1e-12)
    def test_angle(self):
        self.assertAlmostEqual(angle_deg([1,0],[0,1]),90.,places=12)

if __name__=='__main__':unittest.main()
