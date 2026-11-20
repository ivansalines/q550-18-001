"""End-to-end crossing, reconstructed with synthetic 11D sources only."""
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from sign_crossing_math import cumulative_simpson
from run_sign_crossing_v18_7 import bracket_and_rates

class TestSyntheticRunner(unittest.TestCase):
 def test_crossover_from_birth_not_metric(self):
    x=np.linspace(0,8,641);g=np.zeros((len(x),11))
    # Early source points in e0; later births oppose e0 while growing in e1.
    fade=1/(1+np.exp(-40*(x-1.)))
    g[:,0]=1*(1-fade)-.2*fade
    g[:,1]=.12*fade
    I=cumulative_simpson(x,g,1.)
    src={'source_scale':np.array([1.]),'source_s':np.array([0.]),'source_integrand_11D':g}
    mat=lambda s:np.eye(11)
    jac=lambda s,y:np.zeros((11,11))
    state=lambda s:np.zeros(11)
    early=I[80]
    roots=bracket_and_rates('synthetic',src,x,mat,state,jac,I,early)
    self.assertEqual(len(roots),1)
    self.assertGreater(roots[0]['xi_root'],4.)
    self.assertLess(roots[0]['xi_root'],8.)
    self.assertLess(abs(roots[0]['dot_rate_metric']),1e-12)
    self.assertLess(roots[0]['dot_rate_source'],0.)
    self.assertLess(abs(roots[0]['old_total_angle_deg']-90.),1e-7)

if __name__=='__main__':unittest.main()
