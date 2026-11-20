import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from history_window_math import state,orth,eval_at,event_comparison

class TestHistoryAblation(unittest.TestCase):
    def test_constant_source_h_boundary_cancel(self):
        x=np.linspace(0,8,641);e=np.array([1.,0,0])
        A=np.stack([np.zeros_like(x),x,2*x],axis=1)
        g=np.broadcast_to(np.array([0,1,2]),A.shape)
        q=eval_at(x,A,g,1,e,4,.5)
        self.assertAlmostEqual(q['rho'],.5*np.sqrt(5),places=12)
        self.assertAlmostEqual(q['speed'],0,places=12)
        self.assertAlmostEqual(q['radial'],0,places=12)

    def test_moving_window_boundary_is_essential(self):
        x=np.linspace(0,8,641);e=np.array([1.,0,0])
        A=np.stack([np.zeros_like(x),x*x,1+x*x*x],axis=1)
        g=np.stack([np.zeros_like(x),2*x,3*x*x],axis=1)
        z=4.;h=.5
        q=eval_at(x,A,g,1,e,z,h)
        eq=np.array([0.,2*z-2*(z-h),3*z*z-3*(z-h)**2])
        a=np.array([0,z*z-(z-h)**2,z**3-(z-h)**3])
        want=state(e,a,eq)
        self.assertAlmostEqual(q['radial'],want['radial'],places=9)
        self.assertAlmostEqual(q['speed'],want['speed'],places=9)

    def test_old_and_recent_can_be_independent(self):
        x=np.linspace(0,8,641);e=np.array([1.,0,0]);
        A=np.stack([np.zeros_like(x),x,x*x],axis=1)
        g=np.stack([np.zeros_like(x),np.ones_like(x),2*x],axis=1)
        r=event_comparison(x,A,g,1,e,4,.5)
        self.assertEqual(r['rank3'],3)
        self.assertGreater(r['old_exclusive_over_full_radius'],0)

    def test_simple_coordinate_rank_invariant(self):
        a=np.array([1.,0,0,0]);b=np.array([0,1.,1.,0]);c=np.array([0.,0,1.,0]);
        T=np.array([[1.,.3,.2,.1],[0.,2.,.1,0.],[0.,0.,.7,.2],[0.,0.,0.,1.3]])
        self.assertEqual(np.linalg.matrix_rank(np.stack([a,b,c],axis=1)),3)
        self.assertEqual(np.linalg.matrix_rank(np.stack([T@a,T@b,T@c],axis=1)),3)

    def test_degenerate_history_rejected(self):
        with self.assertRaises(ValueError):state(np.array([1.,0.,0.]),np.array([5.,0.,0.]),np.array([0.,1.,0.]))

    def test_lines_do_not_depend_on_sign(self):
        a=state(np.array([1.,0,0]),np.array([0,2.,3]),np.array([0,1,0]))
        b=state(np.array([1.,0,0]),np.array([0,-2.,-3]),np.array([0,-1,0]))
        self.assertAlmostEqual(np.linalg.norm(a['projector']-b['projector']),0,places=12)
        self.assertAlmostEqual(np.linalg.norm(a['pdot']-b['pdot']),0,places=12)

if __name__=='__main__':unittest.main()
