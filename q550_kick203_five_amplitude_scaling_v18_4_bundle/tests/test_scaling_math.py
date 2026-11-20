"""Synthetic, non-Q550 checks of algebra and scaling math only."""
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference'),str(ROOT/'scripts')]
from scaling_math import AMPS,EPS_UNIT,fit_odd
from bridge_math import factor_secants
from run_scan_v18_4 import read_bridge,verify_source

class TestOddFit(unittest.TestCase):
    def test_five_amplitude_configuration(self):
        self.assertTrue(np.array_equal(AMPS,np.array([5e-9,1e-8,1.5e-8,2e-8,3e-8])))

    def test_11d_odd_cubic_quintic(self):
        t=AMPS/EPS_UNIT
        coeff=np.arange(11,dtype=float)+1
        y=t[:,None]*coeff[None,:]+.03*t[:,None]**3*coeff[None,:]-.002*t[:,None]**5*coeff[None,:]
        fit=fit_odd(AMPS,y,powers=(1,3,5))
        self.assertLess(float(np.max(np.abs(fit['residual']))),1e-11)
        self.assertLess(float(np.max(np.abs(fit['loo_residual']))),1e-11)
        self.assertLess(float(np.max(np.abs(fit['coefficients'][0]-coeff))),1e-11)

    def test_cross_factor_cubic_secants(self):
        # 2-factor mixed term is even in the (+/-) central difference and
        # cancels; triple product retains odd cubic cross-factor term.
        eps=1e-3
        descriptor=lambda a,b,c: np.array([a*b*c]+[0.]*10)
        minus=(-eps,-eps,-eps);plus=(eps,eps,eps)
        sec,inter,total=factor_secants(minus,plus,descriptor)
        self.assertLess(np.linalg.norm(sec),1e-20)
        self.assertAlmostEqual(float(inter[0]),2*eps**3,delta=1e-19)
        self.assertLess(np.linalg.norm(total-inter),1e-20)

    def test_v18_3_anchor_data(self):
        d=read_bridge(verify_source())
        self.assertEqual(len(d['xi']),47)
        for tol in ('standard','tight'):
            for eps in (1e-8,3e-8):
                p=f'{tol}__eps{eps:g}'
                delta=d[p+'__centered_delta_11D']
                interaction=d[p+'__centered_interaction_11D']
                sec=d[p+'__centered_factor_secants_3x11']
                self.assertLess(float(np.max(np.linalg.norm(delta-(interaction+sec.sum(axis=1)),axis=1))),1e-12)

if __name__=='__main__':unittest.main()
