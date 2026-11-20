"""Detect deliberately introduced structural corruption before trusting real zeros."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from independent_frozen_audit import block_defects
class Countercheck(unittest.TestCase):
    def test_synthetic_clean(self):
        J=np.zeros((3,11,11));P=np.repeat(np.eye(11)[None],3,axis=0)
        d=block_defects(J,P)
        self.assertEqual(max(d['jac_last_row_max'],d['cocycle_last_row_first10_max'],d['cocycle_last_corner_minus_one_max']),0)
    def test_deliberate_last_row_j_error_is_detected(self):
        J=np.zeros((3,11,11));P=np.repeat(np.eye(11)[None],3,axis=0)
        J[1,10,7]=1e-6
        self.assertGreater(block_defects(J,P)['jac_last_row_max'],1e-9)
    def test_deliberate_last_row_p_error_is_detected(self):
        J=np.zeros((3,11,11));P=np.repeat(np.eye(11)[None],3,axis=0)
        P[1,10,4]=1e-6;P[2,10,10]=1.000001
        d=block_defects(J,P)
        self.assertGreater(d['cocycle_last_row_first10_max'],1e-9)
        self.assertGreater(d['cocycle_last_corner_minus_one_max'],1e-9)
if __name__=='__main__':unittest.main()
