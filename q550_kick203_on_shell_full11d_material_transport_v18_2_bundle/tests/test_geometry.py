#!/usr/bin/env python3
"""Local synthetic controls (NOT original Q550 replay)."""
import sys,unittest
from pathlib import Path
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from material_geometry import (unit,frame_from_flow,nearest_plane_rotation,
                               build_moving_transport,split_increment,angle_deg)

class MaterialGeometryTests(unittest.TestCase):
    def setUp(self):
        self.A=np.zeros((11,11));self.A[0,1]=-1;self.A[1,0]=1
        self.A[1,2]=-.3;self.A[2,1]=.3
        self.initial=np.zeros(11);self.initial[0]=1;self.initial[2]=.25
        self.state=lambda s:expm(s*self.A)@self.initial
        self.rhs=lambda y:self.A@y
        self.jac=lambda y:self.A
    def test_frame_curvature(self):
        E,k=frame_from_flow(self.state(.3),self.rhs,self.jac)
        self.assertEqual(E.shape,(11,2))
        self.assertGreater(k,0)
        self.assertLess(np.linalg.norm(E.T@E-np.eye(2)),1e-10)
        eps=1e-5
        tplus=unit(self.rhs(self.state(.3+eps)))
        tminus=unit(self.rhs(self.state(.3-eps)))
        self.assertLess(np.linalg.norm((tplus-tminus)/(2*eps)-k*E[:,1]),1e-8)
    def test_rotor(self):
        E0,_=frame_from_flow(self.state(0),self.rhs,self.jac)
        E1,_=frame_from_flow(self.state(.2),self.rhs,self.jac)
        for r in (1,2):
            R,a=nearest_plane_rotation(E0[:,:r],E1[:,:r])
            self.assertLess(np.linalg.norm(R@E0[:,:r]@E0[:,:r].T@R.T-E1[:,:r]@E1[:,:r].T),1e-8)
            self.assertEqual(len(a),r)
    def test_transport_and_identity(self):
        samples=[0.,.01,.05,.20,.6]
        for rank in (1,2):
            Rs,diag=build_moving_transport(samples,self.state,self.rhs,self.jac,.045,rank=rank)
            self.assertEqual(len(Rs),len(samples))
            self.assertGreater(diag['steps'],len(samples)-1)
            u0=self.state(.05)+np.ones(11)*.03
            u1=self.state(.2)+np.ones(11)*.02
            internal,frame,error=split_increment(u0,u1,Rs[.05],Rs[.2])
            self.assertLess(error,1e-10)
            self.assertLess(np.linalg.norm(internal+frame-(u1-u0)),1e-10)
    def test_no_physical_flip_from_rigid_turn(self):
        # A single moving vector can rotate in the fixed lab coordinates while
        # its one-dimensional co-transported direction remains unchanged.
        E0,_=frame_from_flow(self.state(0),self.rhs,self.jac)
        E1,_=frame_from_flow(self.state(.05),self.rhs,self.jac)
        R,_=nearest_plane_rotation(E0[:,:1],E1[:,:1])
        u0=E0[:,0];u1=R@u0
        self.assertLess(np.linalg.norm(u1-R@u0),1e-11)
        self.assertAlmostEqual(float((R.T@u1)@u0),1.,places=9)
    def test_degeneracy_guard(self):
        with self.assertRaises(ValueError):
            frame_from_flow(np.zeros(11),self.rhs,self.jac)
        self.assertTrue(np.isnan(angle_deg(np.zeros(11),np.ones(11))))

if __name__=='__main__':unittest.main()
