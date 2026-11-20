#!/usr/bin/env python3
import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from reorientation_math import quantities,sphere_tangent_parallel_transport,coherence,roots,root_pair_coherence

class TestGeometry(unittest.TestCase):
    def test_radial_event_continued_rotation(self):
        x=np.linspace(-.5,.5,101);rho=2+x*x
        e=np.array([1.,0,0]);v=np.stack([np.zeros_like(x),np.cos(x),np.sin(x)],axis=1)
        A=v*rho[:,None]
        g=np.gradient(A,x,axis=0,edge_order=2)
        q=quantities(e,A,g,1)
        self.assertTrue(np.any(q['radial_rate'][:-1]*q['radial_rate'][1:]<=0))
        self.assertAlmostEqual(float(q['speed'][50]),1,places=3)
        p=root_pair_coherence(x,e,A,g,1,0,0.10)
        self.assertGreater(p['cosine'],.9999)

    def test_tangent_flip_control(self):
        # f(x)=x^2 generates opposite tangent directions around zero, unlike a smoothly
        # continuing angular passage; the detector must report anti-coherence.
        x=np.linspace(-.5,.5,1001);rho=2+x*x;ang=x*x
        e=np.array([1.,0,0]);v=np.stack([np.zeros_like(x),np.cos(ang),np.sin(ang)],axis=1)
        A=v*rho[:,None];g=np.gradient(A,x,axis=0,edge_order=2)
        p=root_pair_coherence(x,e,A,g,1,0,0.1)
        self.assertLess(p['cosine'],-.9999)

    def test_spherical_parallel_transport_orthogonality(self):
        a=np.array([1.,0,0]);b=np.array([np.cos(.7),np.sin(.7),0.]);t=np.array([0.,1,0])
        moved=sphere_tangent_parallel_transport(a,b,t)
        self.assertAlmostEqual(float(moved@b),0,places=12)
        self.assertAlmostEqual(float(np.linalg.norm(moved)),1,places=12)
        self.assertGreater(coherence(a,t,b,moved),.9999999)

    def test_zero_source_no_rotation(self):
        e=np.array([1.,0.,0.]);a=np.array([[0.,3,0],[0,3,0]]);g=np.zeros((2,3))
        self.assertAlmostEqual(float(quantities(e,a,g,1)['speed'].sum()),0.)

    def test_metric_scaling_invariance(self):
        e=np.array([1.,2.,.5]);a=np.array([[1.,3.,2.]]);g=np.array([[3.,-1.,2.]])
        q1=quantities(e,a,g,1);q2=quantities(e*2,a*2,g*2,1)
        self.assertAlmostEqual(float(q1['speed'][0]),float(q2['speed'][0]),places=13)
        self.assertAlmostEqual(float(q1['fraction_outside'][0]),float(q2['fraction_outside'][0]),places=13)

    def test_root_interpolation(self):
        x=np.linspace(-.4,.4,20);v=x-0.123456
        for mode in ['linear','pchip','cubic']:
            self.assertAlmostEqual(roots(x,v,mode)[0],.123456,places=10)

    def test_nonzero_wedge_rank_not_metric(self):
        e=np.array([1.,0,0,0]);a=np.array([0.,1,0,0]);g=np.array([0.,0,1,0]);
        T=np.array([[2.,.4,.1,0],[0.,.4,0.,1.],[0.,0.,3,0.],[0.,0,0,1.]])
        self.assertEqual(np.linalg.matrix_rank(np.stack([e,a,g],axis=1)),3)
        self.assertEqual(np.linalg.matrix_rank(np.stack([T@e,T@a,T@g],axis=1)),3)

if __name__=='__main__':unittest.main()
