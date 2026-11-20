#!/usr/bin/env python3
"""Analytic, no-server tests: second-order differences and frozen genealogy."""
import io
import sys
import unittest
import zipfile
from pathlib import Path

import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'reference')]
from run_kick203_second_order_observer_v18_1 import tangent_differences,read_npz_from_zip,angle_deg
from frame_geometry import shapley_three_factors

class SmoothHistory:
    def tangent(self,s):
        v=np.array([np.cos(s),np.sin(s),0.2*np.sin(2*s),0.15*np.cos(3*s)]+[0.] *7)
        return v / np.linalg.norm(v)

class Tests(unittest.TestCase):
    def test_frame_orthonormal(self):
        for h in [.02,.01,.005,.0025]:
            frames=tangent_differences(SmoothHistory(),.61,h)
            for f in frames.values():
                self.assertLess(np.linalg.norm(f.T@f-np.eye(2)),1e-12)

    def test_unilateral_second_order_converges_faster(self):
        tiny=tangent_differences(SmoothHistory(),.61,1e-5)['centered']
        for direction in ('forward','backward'):
            f1=tangent_differences(SmoothHistory(),.61,.02)
            f2=tangent_differences(SmoothHistory(),.61,.01)
            e1=np.linalg.norm(f1[direction+'1']-tiny)
            e1half=np.linalg.norm(f2[direction+'1']-tiny)
            e2=np.linalg.norm(f1[direction+'2']-tiny)
            e2half=np.linalg.norm(f2[direction+'2']-tiny)
            self.assertGreater(e1 / e1half,1.8)
            self.assertGreater(e2 / e2half,3.2)
            self.assertLess(e2,e1)

    def test_shapley_identity_11d(self):
        rng=np.random.default_rng(12)
        a=np.linalg.qr(rng.normal(size=(11,2)))[0]
        b=np.linalg.qr(rng.normal(size=(11,2)))[0]
        M=np.eye(11)
        res,delta,err=shapley_three_factors((a,b,M),(a+.001*rng.normal(size=(11,2)),b+.001*rng.normal(size=(11,2)),M+.001*rng.normal(size=(11,11))))
        self.assertLess(err,1e-13)
        self.assertLess(np.linalg.norm(res.sum(axis=0)-delta),1e-13)

    def test_archive_shape_and_keys(self):
        p=ROOT/'inputs/q550_kick203_one_sided_observer_v18_0_sendback.zip'
        d=read_npz_from_zip(p,'observer_shapley_arrays.npz')
        xi=d['xi']
        self.assertEqual(len(xi),47)
        self.assertTrue(np.all(np.diff(xi)>0))
        k='standard__kick203__eps1e-08__hquarter__non_straddling__dlocal_delta'
        self.assertEqual(d[k].shape,(len(xi),11))
        self.assertEqual(np.sum(~np.isfinite(d[k]).all(axis=1)),1)
        self.assertTrue(np.isfinite(d[k][np.isfinite(d[k]).all(axis=1)]).all())

if __name__=='__main__': unittest.main()
