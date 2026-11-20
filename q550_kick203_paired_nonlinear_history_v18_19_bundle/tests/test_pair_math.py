import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from pair_math import chosen_seed,orth_rows
class MathTests(unittest.TestCase):
 def setUp(self):
  self.h=np.eye(10)[:4]; self.j=np.random.default_rng(1819).normal(size=(10,10))
 def test_hidden_jet(self):
  v,d=chosen_seed(self.h,self.j,'jet')
  self.assertLess(np.linalg.norm(self.h@v),1e-10)
  self.assertLess(np.linalg.norm(self.h@self.j@v),1e-10)
 def test_hidden_state_only(self):
  v,d=chosen_seed(self.h,self.j,'state_only')
  self.assertLess(np.linalg.norm(self.h@v),1e-10)
  self.assertGreater(np.linalg.norm(self.h@self.j@v),.1)
 def test_observer_rotation_covariance(self):
  A=np.random.default_rng(2026).normal(size=(10,3));h=orth_rows(A)
  self.assertLess(np.linalg.norm(h@h.T-np.eye(3)),1e-12)
 def test_both_signs_equal_reduced_state(self):
  v,d=chosen_seed(self.h,self.j,'jet')
  p=np.ones(10)+1e-4*v;m=np.ones(10)-1e-4*v
  self.assertLess(np.linalg.norm(self.h@(p-m)),1e-10)
 def test_matched_derivative_is_not_proof_of_identity(self):
  v,_=chosen_seed(self.h,self.j,'jet')
  self.assertGreater(np.linalg.norm(v),.99)
  self.assertLess(np.linalg.norm(self.h@v)+np.linalg.norm(self.h@self.j@v),1e-10)
 def test_source_model_not_certified_as_physical(self):
  # Test insists no physical property is inferred from linear algebra.
  v,d=chosen_seed(self.h,self.j,'jet');self.assertGreater(d['hidden_dim'],0)
if __name__=='__main__':unittest.main()

class SyntheticHistoryTest(unittest.TestCase):
 def test_equal_observable_and_jet_can_diverge_later(self):
  # x'=y, y'=z, z'=0; x(0) and x'(0) equal for z=+/-eps.
  # x_+(t)-x_-(t)=eps*t^2, arising from a HIDDEN PRESENT state,
  # not constitutive memory. Independent nonlinear integrations.
  from scipy.integrate import solve_ivp
  def rhs(t,y): return [y[1],y[2],0.]
  e=2e-5
  solp=solve_ivp(rhs,(0,2),[0.,0.,e],rtol=1e-12,atol=1e-14,t_eval=[0.,1.,2.])
  solm=solve_ivp(rhs,(0,2),[0.,0.,-e],rtol=1e-12,atol=1e-14,t_eval=[0.,1.,2.])
  self.assertEqual(solp.y[0,0],solm.y[0,0])
  self.assertEqual(solp.y[1,0],solm.y[1,0])
  self.assertAlmostEqual(solp.y[0,-1]-solm.y[0,-1],e*4,places=12)
 def test_equal_full_state_equal_future(self):
  from scipy.integrate import solve_ivp
  fun=lambda t,y:[y[1],-.2*y[0]]
  args=((0.,1.),[.3,-.2])
  a=solve_ivp(fun,*args,t_eval=[.25,.5,1.],rtol=1e-12,atol=1e-14)
  b=solve_ivp(fun,*args,t_eval=[.25,.5,1.],rtol=1e-12,atol=1e-14)
  self.assertLess(np.max(np.abs(a.y-b.y)),1e-13)
