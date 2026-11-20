#!/usr/bin/env python3
"""No-server unit tests for central/forward/backward difference geometry."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path[:0]=[str(Path(__file__).resolve().parents[1]/'reference'),str(Path(__file__).resolve().parents[1]/'scripts')]
from run_kick203_one_sided_observer_v18_0 import three_frames
from frame_geometry import shapley_three_factors
class TestOneSidedFrame(unittest.TestCase):
 def test_agree_for_linear_unit_tangent(self):
  class Fake:
   def tangent(self,s):
    t=np.array([np.cos(0.001*s),np.sin(0.001*s)]+[0.]*9)
    return t
  d=three_frames(Fake(),1.0,.01)
  self.assertLess(np.linalg.norm(d['centered']-d['forward']),.00001)
  self.assertLess(np.linalg.norm(d['centered']-d['backward']),.00001)
 def test_transverse_nonlinear_shift(self):
  class Fake:
   def tangent(self,s):
    a=np.array([1.,s,s*s]+[0.]*8)
    return a/np.linalg.norm(a)
  d=three_frames(Fake(),.5,.1)
  self.assertGreater(np.linalg.norm(d['forward']-d['backward']),.05)
  for f in d.values():self.assertLess(np.linalg.norm(f.T@f-np.eye(2)),1e-12)
 def test_shapley_11d_closure(self):
  rng=np.random.default_rng(42)
  A=np.eye(11)[:,:2];B=np.eye(11)[:,2:4];M=np.eye(11)
  # A and B must not be orthogonal for meaningful descriptor variation.
  B=(A*.5 +B*.8660254038)
  a=(A,B,M);b=(A+.001*rng.normal(size=(11,2)),B+.002*rng.normal(size=(11,2)),M+.001*rng.normal(size=(11,11)))
  share,delta,err=shapley_three_factors(a,b)
  self.assertLess(err,1e-12)
  self.assertTrue(np.allclose(share.sum(axis=0),delta))
if __name__=='__main__':unittest.main()
