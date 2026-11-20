"""Synthetic tests for the v18.3 finite-history identities and input gates."""
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'reference'))
from bridge_math import factor_secants, match_grid
from history_bridge_v16 import descriptor,gram_frame
from frame_geometry import shapley_three_factors

class TestBridge(unittest.TestCase):
    def test_conditional_exact_nonlinear(self):
        def d(a,b,c):return np.array([a[0]*b[0]+np.sin(c[0]),a[0]*b[0]*c[0]])
        m=(np.array([1.]),np.array([2.]),np.array([3.]))
        p=(np.array([4.]),np.array([3.]),np.array([2.]))
        sec,inter,total=factor_secants(m,p,d)
        self.assertTrue(np.allclose(np.sum(sec,axis=0)+inter,total,atol=1e-15))
        self.assertGreater(np.linalg.norm(inter),0.0)
    def test_exact_frame_descriptor(self):
        rng=np.random.default_rng(550)
        X=rng.normal(size=(11,11))
        q,_=np.linalg.qr(X)
        a0=gram_frame(q[:,0],q[:,1]);b0=gram_frame(q[:,2],q[:,3])
        m0=np.eye(11)
        a1=gram_frame(q[:,0]+.03*q[:,4],q[:,1]+.02*q[:,5])
        b1=gram_frame(q[:,2]+.05*q[:,6],q[:,3]-.01*q[:,8])
        m1=np.eye(11)+.01*np.outer(q[:,4],q[:,1])
        lo=(a0,b0,m0);hi=(a1,b1,m1)
        sec,interaction,delta=factor_secants(lo,hi,descriptor)
        share,exact,closure=shapley_three_factors(lo,hi)
        self.assertLess(np.linalg.norm(sum(sec)+interaction-delta),1e-13)
        self.assertLess(np.linalg.norm(sum(share)-delta),1e-13)
        self.assertLess(closure,1e-13)
        self.assertLess(np.linalg.norm(delta-exact),1e-13)
    def test_alignment(self):
        a=np.array([0.,.01,.02,.03]);b=np.array([.01,.03])
        self.assertTrue(np.array_equal(match_grid(a,b),[1,3]))
        with self.assertRaises(ValueError):match_grid(a,[.012])
        with self.assertRaises(ValueError):match_grid(a,[.03,.01])
    def test_factor_rejection(self):
        with self.assertRaises(ValueError):factor_secants((np.ones(1),np.ones(2),np.ones(3)),(np.ones(3),np.ones(2),np.ones(3)),lambda *x: x[0])

if __name__=='__main__':unittest.main()
