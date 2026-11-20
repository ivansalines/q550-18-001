"""End-to-end figure/report postprocessing with synthetic Q550-shaped arrays."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import run_sign_crossing_v18_7 as runner
import analyze_sign_crossing_v18_7 as report
from sign_crossing_math import cumulative_simpson

class TestReporting(unittest.TestCase):
 def test_reports_and_plots(self):
  xi=np.linspace(0.,8.,641);fade=1/(1+np.exp(-40*(xi-1)))
  g=np.zeros((641,11));g[:,0]=1-1.2*fade;g[:,1]=.12*fade
  src={'source_s':np.array([0]),'source_scale':np.array([1]),'xi_nodes':xi,
       'source_integrand_11D':g,'source_integrand_halfstep_11D':g,
       'cocycle_condition_nodes':np.ones(641),
       'predicted_quadratic_at_delays_11D':cumulative_simpson(xi,g,1.)[[20,40,80,160,320,640]]}
  class P:
   def sol(self,t):return np.eye(11).ravel()
  with TemporaryDirectory() as tmp:
   saved=(runner.ROOT,report.ROOT,report.OUT)
   try:
    runner.ROOT=Path(tmp);report.ROOT=Path(tmp);report.OUT=Path(tmp)/'results'
    for tol in ('standard','tight'):
     runner.analyze_tol(tol,src,P(),lambda s:np.zeros(11),
           lambda s,y:np.zeros((11,11)),SimpleNamespace(condition_gate=1e8,replay_relative_gate=1e-8))
    report.main()
    for file in ['REPORT.md','summary.json','source_birth_sign_crossing_dense.png',
                 'crossing_rate_source_vs_transport.png','polar_corotation_and_source_frame_angles.png']:
     self.assertTrue((Path(tmp)/'results'/file).is_file(),file)
   finally:runner.ROOT,report.ROOT,report.OUT=saved
if __name__=='__main__':unittest.main()
