#!/usr/bin/env python3
"""Synthetic algebra-only tests. Synthetic outcomes are never Q550 measurements."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from structural_math_v18_9 import (roots_bracketed,cosines,transport_invariant,
    symmetry_basis,lyapunov_null_screen,affine_invariant_screen,unit_change_check,euler_strain_birth)

class StructuralTests(unittest.TestCase):
    def test_double_sign_change(self):
        x=np.linspace(0,8,81); y=(x-2.03)*(x-5.07); yp=2*x-7.10
        got=roots_bracketed(x,y,yp)
        self.assertEqual(len(got),2)
        self.assertLess(max(abs(got-np.array([2.03,5.07]))),1e-10)
    def test_source_frame_transport(self):
        rng=np.random.default_rng(9);P=np.stack([np.diag(np.linspace(0.8,1.7,11)) for _ in range(5)])
        E=rng.normal(size=11);T=rng.normal(size=(5,11))
        got=transport_invariant(P,E,T)
        self.assertTrue(np.allclose(got,T@E,atol=1e-12))
    def test_tensor_covariance_under_units(self):
        rng=np.random.default_rng(9);E=rng.normal(size=(5,11));K=rng.normal(size=(5,11));P=np.broadcast_to(np.eye(11),(5,11,11))
        self.assertLess(unit_change_check(P,E,K,np.exp(np.arange(11)/5),np.ones(11)),1e-12)
    def test_symmetric_screen_has_66_directions(self):
        B,pairs=symmetry_basis();self.assertEqual(B.shape,(66,11,11))
        J=np.zeros((16,11,11));J[:,0,0]=-1
        screen=lyapunov_null_screen(J,sample_stride=2)
        self.assertEqual(screen['n_symmetric_unknowns'],66)
        screen2=affine_invariant_screen(J,sample_stride=2)
        self.assertAlmostEqual(screen2['singular_min'],0)
    def test_strain_and_source_split(self):
        J=np.zeros((3,11,11));P=np.stack([np.eye(11)]*3);E=np.ones((3,11));K=E.copy();g=K.copy()
        st,so=euler_strain_birth(np.arange(3),E,K,J,P,g,.1)
        self.assertTrue(np.allclose(st,0));self.assertTrue(np.allclose(so,1.1))
    def test_metric_can_change_sign_independent_vectors(self):
        e=np.array([1.,0.]);k=np.array([0.,1.]);
        self.assertLess(e@np.array([[2.,-1.],[-1.,2.]])@k,0)
        self.assertGreater(e@np.array([[2.,1.],[1.,2.]])@k,0)

if __name__=='__main__':unittest.main()
