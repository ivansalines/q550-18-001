import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from residual_plane_math import source_geometry,principal_line_angle, max_tensor_transform_error

class Tests(unittest.TestCase):
    def test_pure_transport_source_zero(self):
        e=np.array([1.,0,0]);a=np.array([[0.,2.,0.],[0.,2.,0.]]);g=np.zeros((2,3))
        r=source_geometry(e,a,g,2)
        self.assertTrue(np.allclose(r['turn_rate'],0))
    def test_source_parallel_to_current_plane_cannot_turn(self):
        e=np.array([1.,0,0]);a=np.array([[0.,1.,0.]]);g=np.array([[3.,2.,0.]])
        self.assertAlmostEqual(source_geometry(e,a,g,1)['turn_rate'][0],0)
    def test_outside_source_must_turn(self):
        e=np.array([1.,0,0]);a=np.array([[0.,1.,0.]]);g=np.array([[0.,0.,4.]])
        self.assertAlmostEqual(source_geometry(e,a,g,.5)['turn_rate'][0],2)
    def test_nonuniform_coordinate_change_raw_wedge(self):
        e=np.array([1.,2,3]);a=np.array([3.,1,2]);g=np.array([3.,4,5]);T=np.array([[3.,1.,0],[0.,2.,1.],[0,0,.5]])
        self.assertLess(max_tensor_transform_error(T,e,a,g),1e-13)
    def test_line_not_vector_angle(self):
        v=np.array([[1.,0,0.]]);w=np.array([[-1.,0,0.]])
        self.assertAlmostEqual(float(principal_line_angle(v,w)[0]),0)
    def test_degeneracy_detected(self):
        with self.assertRaises(ValueError):source_geometry(np.array([1.,0.,0.]),np.array([[3.,0,0.]]),np.array([[0,1.,0]]),1)

if __name__=='__main__':unittest.main()
