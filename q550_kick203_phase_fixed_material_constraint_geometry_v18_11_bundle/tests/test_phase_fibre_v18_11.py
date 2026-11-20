"""Numerical synthetic tests for phase block-cocycle and Grassmann derivative."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from phase_fibre_math_v18_11 import phase_sector_checks,phase_cocycle_rhs,projector_derivative,principal_step_angle
class TestPhaseFibre(unittest.TestCase):
 def test_invariant_vertical_subbundle(self):
  J=np.zeros((3,11,11));P=np.tile(np.eye(11),(3,1,1));J[:,:10,10]=2
  record=phase_sector_checks(J,P,[np.zeros((3,11))])
  self.assertEqual(record['max_jacobian_phase_row'],0)
  self.assertEqual(record['max_cocycle_phase_row_to_z'],0)
  self.assertEqual(record['max_phase_to_z_jacobian_column_norm'],2*np.sqrt(10))
 def test_forced_phase_variational_block(self):
  J=np.zeros((2,11,11));P=np.tile(np.eye(11),(2,1,1));J[:,0,0]=3;J[:,1,10]=4
  Phi,q=phase_cocycle_rhs(J,P,2)
  self.assertAlmostEqual(Phi[0,0,0],6);self.assertAlmostEqual(q[0,1],8)
 def test_projector_rate_against_difference(self):
  rng=np.random.default_rng(1811);U=rng.normal(size=(10,2));dU=rng.normal(size=(10,2));h=1e-6
  Q,_,_=np.linalg.svd(U,full_matrices=False); Q=Q[:,:2]
  fun=lambda Z: (lambda W:W@W.T)(np.linalg.svd(Z,full_matrices=False)[0])
  numerical=(fun(U+h*dU)-fun(U-h*dU))/(2*h)
  derivative,*_=projector_derivative(U,dU)
  np.testing.assert_allclose(numerical,derivative,atol=3e-9,rtol=3e-9)
 def test_line_plane_angle(self):
  U=np.eye(10)[:,:2]
  self.assertAlmostEqual(principal_step_angle(U,U),0.,places=6)
  V=U.copy();theta=.03;V[0,0]=np.cos(theta);V[2,0]=np.sin(theta)
  self.assertAlmostEqual(principal_step_angle(U,V),theta,places=6)
 def test_source_birth_additivity(self):
  rng=np.random.default_rng(11);U=rng.normal(size=(10,2));dt=rng.normal(size=(10,2));db=rng.normal(size=(10,2))
  all_,*_=projector_derivative(U,dt+db);t,*_=projector_derivative(U,dt);b,*_=projector_derivative(U,db)
  np.testing.assert_allclose(all_,t+b,atol=1e-12,rtol=1e-12)
if __name__=='__main__':unittest.main()
