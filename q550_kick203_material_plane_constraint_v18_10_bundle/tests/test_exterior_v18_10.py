"""Exact exterior algebra synthetic invariance and 11D transport tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from exterior_geometry_v18_10 import blade2,blade3,wedge_rate,rank_from_vectors
class ExteriorTests(unittest.TestCase):
 def setUp(self):
  self.rng=np.random.default_rng(1810)
 def test_antisymmetry_and_collinearity(self):
  a=self.rng.normal(size=11);b=self.rng.normal(size=11)
  self.assertLess(np.linalg.norm(blade2(a,a)),1e-12)
  np.testing.assert_allclose(blade2(a,b),-blade2(b,a),rtol=1e-13,atol=1e-13)
  self.assertLess(np.linalg.norm(blade3(a,b,a)),1e-12)
 def test_algebraic_invertible_covariance(self):
  a,b,c=self.rng.normal(size=(3,11));T=self.rng.normal(size=(11,11));T+=5*np.eye(11)
  self.assertEqual(rank_from_vectors([a,b,c])[0],rank_from_vectors([T@a,T@b,T@c])[0])
  self.assertEqual(rank_from_vectors([a,b,a+b])[0],rank_from_vectors([T@a,T@b,T@(a+b)])[0])
  self.assertLess(np.linalg.norm(blade3(T@a,T@b,T@(a+b))),1e-11)
 def test_continuous_transport_rate(self):
  J=np.diag(np.linspace(-.1,.2,11));u=self.rng.normal(size=11);a=self.rng.normal(size=11)
  w=self.rng.normal(size=11);xi=.43;g=w
  def obs(x):
   P=np.diag(np.exp(x*np.diag(J)));A=a+x*w
   return P@u,P@A,P
  E,K,P=obs(xi); h=1e-6
  ep,kp,_=obs(xi+h);em,km,_=obs(xi-h)
  numeric=(blade2(ep,kp)-blade2(em,km))/(2*h)
  strain,birth=wedge_rate(E[None,:],K[None,:],J[None,:,:],P[None,:,:],g[None,:],1.)
  np.testing.assert_allclose(numeric,strain[0]+birth[0],atol=1e-9,rtol=1e-8)
 def test_bivector_is_not_metric_dot(self):
  a=np.zeros(11);b=np.zeros(11);a[0]=1;b[1]=1
  self.assertAlmostEqual(float(a@b),0.)
  self.assertGreater(np.linalg.norm(blade2(a,b)),.9)
 def test_rank3_outside_plane(self):
  a=np.eye(11)[0];b=np.eye(11)[1];g=np.eye(11)[2]
  self.assertEqual(rank_from_vectors([a,b,g])[0],3)
  self.assertEqual(rank_from_vectors([a,b,a+3*b])[0],2)
if __name__=='__main__':unittest.main()
