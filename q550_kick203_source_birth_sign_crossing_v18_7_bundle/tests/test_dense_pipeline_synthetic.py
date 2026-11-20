"""Integration test: dense source data through all analysis and controls (synthetic)."""
import sys,unittest,tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import run_sign_crossing_v18_7 as audit
from sign_crossing_math import cumulative_simpson

class TestSyntheticDensePipeline(unittest.TestCase):
 def test_full_npz_and_crossing_controls(self):
  xi=np.linspace(0,8,641);g=np.zeros((641,11));fade=1/(1+np.exp(-40*(xi-1)))
  g[:,0]=1-fade-.2*fade;g[:,1]=.12*fade
  integ=cumulative_simpson(xi,g,1)
  sample_idx=np.asarray([20,40,80,160,320,640]);targets=integ[sample_idx]
  src={'source_s':np.array([0.]),'source_scale':np.array([1.]),
       'xi_nodes':xi,'source_integrand_11D':g,'source_integrand_halfstep_11D':g,
       'predicted_quadratic_at_delays_11D':targets,'cocycle_condition_nodes':np.ones(641)}
  class Fundamental:
   def sol(self,t):return np.eye(11).ravel()
  with tempfile.TemporaryDirectory() as tmp:
   save=audit.ROOT
   try:
    audit.ROOT=Path(tmp)
    result=audit.analyze_tol('synthetic',src,Fundamental(),lambda s:np.zeros(11),
       lambda s,y:np.zeros((11,11)),SimpleNamespace(condition_gate=1e8,replay_relative_gate=1e-8))
    self.assertEqual(len(result['roots']),1)
    self.assertLess(abs(result['roots'][0]['dot_rate_metric']),1e-12)
    self.assertLess(result['roots'][0]['dot_rate_source'],0)
    self.assertTrue((Path(tmp)/'checkpoint/synthetic_sign_crossing.npz').exists())
    self.assertEqual(len((Path(tmp)/'results/synthetic_dense_projection.csv').read_text().splitlines()),562)
   finally:audit.ROOT=save

if __name__=='__main__':unittest.main()
