import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from jet_math import derivatives,motion_matrices,observability,numerical_rank,null_vectors,canonical,scaled_taylor
class JetMathTests(unittest.TestCase):
 def test_constant_derivatives(self):
  a=np.array([[1.,2.],[-3.,4.]])
  aa,d1,d2=derivatives(np.repeat(a[None,:,:],25,axis=0),12,.1)
  self.assertLess(np.max(abs(d1)),1e-12)
  self.assertLess(np.max(abs(d2)),1e-12)
  m=motion_matrices(aa,d1,d2)
  self.assertLess(np.max(abs(m[3]-a@a@a)),1e-12)
 def test_quadratic_time_coefficients(self):
  a=np.array([[.1,1.],[0.,-.5]])
  b=np.array([[.4,0.],[2.,0.]])
  c=np.array([[0.,.3],[.7,0.]])
  h=.0125;xs=np.arange(-20,21)*h
  arr=np.array([a+x*b+x*x*c/2 for x in xs]);aa,d1,d2=derivatives(arr,20,h)
  self.assertLess(np.max(abs(d1-b)),1e-11)
  self.assertLess(np.max(abs(d2-c)),1e-10)
 def test_exact_taylor_scalar(self):
  a0=.4; b=.2; c=-.3;h=.025
  A=np.array([[[a0+x*b+.5*x*x*c]] for x in (np.arange(21)-10)*h])
  aa,d1,d2=derivatives(A,10,h);m=motion_matrices(aa,d1,d2)
  exact=lambda x:np.exp(a0*x+b*x*x/2+c*x*x*x/6)
  # Cubic Taylor residual must scale ~h^4, not exactly zero.
  err=[]
  for x in [.05,.1]:
   approx=scaled_taylor(np.ones((1,1)),m,np.ones(1),x)[0]
   err.append(abs(exact(x)-approx))
  self.assertTrue(10<err[1]/err[0]<25)
 def test_nullspace(self):
  H=np.array([[1.,0.,0.],[0.,1.,0.]])
  N=null_vectors(H);self.assertEqual(N.shape,(3,1))
  self.assertLess(np.linalg.norm(H@N),1e-12)
 def test_rank(self):
  X=np.diag([1.,.05,1e-11]);rk,s=numerical_rank(X)
  self.assertEqual(rk,2)
 def test_equivariant_change_coordinates(self):
  A=np.array([[1.,2.],[-1.,.4]]);b=np.eye(2)*.2;d=np.eye(2)*.5
  T=np.array([[2.,.3],[0.,1.7]])
  m=motion_matrices(A,b,d)
  n=motion_matrices(T@A@np.linalg.inv(T),T@b@np.linalg.inv(T),T@d@np.linalg.inv(T))
  for x,y in zip(m,n):self.assertLess(np.linalg.norm(T@x@np.linalg.inv(T)-y),1e-10)
 def test_jet_rank_example(self):
  A=np.array([[0.,1.,0.],[0.,0.,1.],[0.,0.,0.]])
  M=motion_matrices(A,np.zeros_like(A),np.zeros_like(A))
  O=observability(np.array([[1.,0.,0.]]),M)
  self.assertEqual([numerical_rank(q)[0] for q in O],[1,2,3,3])
 def test_hidden_seed_jet_conditions(self):
  A=np.array([[0.,1.,0.],[0.,0.,1.],[0.,0.,0.]])
  M=motion_matrices(A,np.zeros_like(A),np.zeros_like(A))
  H=np.array([[1.,0.,0.]])
  N=null_vectors(observability(H,M)[1]);v=canonical(N[:,0])
  self.assertLess(np.linalg.norm(H@M[0]@v),1e-12)
  self.assertLess(np.linalg.norm(H@M[1]@v),1e-12)
  self.assertGreater(np.linalg.norm(H@M[2]@v),.9)
 def test_stencil_bounds(self):
  with self.assertRaises(ValueError):derivatives(np.zeros((7,2,2)),1,.1,1)
if __name__=='__main__':unittest.main()
