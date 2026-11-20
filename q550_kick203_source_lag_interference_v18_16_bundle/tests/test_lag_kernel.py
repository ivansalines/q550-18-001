#!/usr/bin/env python3
"""Unit tests: mathematically independent toy sources and local invariants."""
import unittest,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from lag_kernel_math import project,age_slabs,kernel,kernel_from_slabs,pairing_identity
class TestKernel(unittest.TestCase):
    def test_constant_source_has_unit_kernel(self):
        q=np.tile([2.,-3.],(20,1))
        for k in [0,1,5,15]:self.assertAlmostEqual(kernel_from_slabs(q,k,np.ones(2))[0],1.,places=13)
    def test_alternating_source_negative_lag_one(self):
        q=np.array([[(-1)**j,0.] for j in range(20)])
        self.assertAlmostEqual(kernel_from_slabs(q,1,np.ones(2))[0],-1.,places=13)
        self.assertAlmostEqual(kernel_from_slabs(q,2,np.ones(2))[0],1.,places=13)
    def test_metric_project_exact(self):
        e=np.arange(1.,11.);w=np.geomspace(.5,2.,10)
        q=project(e,np.arange(10.,0.,-1.),w)
        self.assertAlmostEqual(float(q@(w*e)),0.,places=12)
    def test_additive_source_closure(self):
        x=np.linspace(0,8,641);A=np.vstack([x,x*x,np.sin(x)]).T
        ai=PchipInterpolator(x,A,axis=0);e=ai(1.)
        q,length=age_slabs(ai,e,4.,.25,np.ones(3))
        self.assertAlmostEqual(length,3.)
        self.assertLess(np.linalg.norm(q.sum(0)-project(e,ai(4.)-ai(1.),np.ones(3))),1e-12)
        self.assertLess(pairing_identity(q,np.ones(3)),1e-10)
    def test_support_rejection(self):
        x=np.linspace(0,8,641);A=np.column_stack([x,np.sin(x)])
        ai=PchipInterpolator(x,A,axis=0);rr,covered=kernel(ai,ai(1.),2.,.125,[.5,1.5],np.ones(2),None,.5)
        self.assertTrue(np.isfinite(rr[0][1]));self.assertTrue(np.isnan(rr[1][1]))
    def test_bounded(self):
        rng=np.random.default_rng(44)
        q=rng.normal(size=(80,10));w=np.linspace(.8,1.2,10)
        for k in range(40):self.assertLessEqual(abs(kernel_from_slabs(q,k,w)[0]),1+1e-12)
    def test_coordinate_rescaling_with_metric(self):
        rng=np.random.default_rng(7);q=rng.normal(size=(40,10));w=np.geomspace(.5,2.,10)
        t=np.geomspace(.6,1.8,10);qp=q*t;wp=w/(t*t)
        for k in (1,8,17):self.assertAlmostEqual(kernel_from_slabs(q,k,w)[0],kernel_from_slabs(qp,k,wp)[0],places=13)
if __name__=='__main__':unittest.main()
