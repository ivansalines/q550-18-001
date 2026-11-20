"""Synthetic invariants and exact algebra tests; NEVER represent Q550 results."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from metric_math_v18_8 import (projectors_from_jacobians,fixed_metric_dot_and_rate,
 source_metric_dot_and_rate,roots_from_samples,mixed_metric,dot_derivative_variable_metric,
 fixed_w_from_reference)

class MetricTests(unittest.TestCase):
    def test_symmetric_projector(self):
        j=np.zeros((5,2,11));j[:,0,0]=1;j[:,1,2]=3
        proj,r,g=projectors_from_jacobians(j)
        self.assertEqual(r,2)
        np.testing.assert_allclose(proj@proj,proj,atol=1e-13)
        np.testing.assert_allclose(np.einsum('nmk,nkj->nmj',j,proj),0.,atol=1e-12)
    def test_rank_variable_rejected(self):
        j=np.zeros((3,1,11));j[0,0,0]=1
        with self.assertRaises(ValueError):projectors_from_jacobians(j)
    def test_unit_transform_invariance_with_congruence(self):
        rng=np.random.default_rng(10);e=rng.normal(size=11);k=rng.normal(size=11)
        W=np.eye(11);D=np.diag(np.linspace(.4,3.2,11));
        self.assertAlmostEqual(e@W@k,(D@e)@np.linalg.inv(D).T@W@np.linalg.inv(D)@(D@k),places=11)
    def test_fixed_rate_finite_difference_on_linear_flow(self):
        rng=np.random.default_rng(8);j=np.diag(np.linspace(-.2,.3,11));g=np.linspace(.1,.3,11)
        x=np.linspace(0,2,201);P=np.array([np.diag(np.exp(np.diag(j)*v)) for v in x]);
        a=np.ones(11);T=np.array([a+.5*v*g for v in x]);
        E=np.einsum('nij,j->ni',P,a);K=np.einsum('nij,nj->ni',P,T)
        gg=np.tile(.5*g,(len(x),1));JJ=np.tile(j,(len(x),1,1)); W=np.diag(np.linspace(.5,2.,11))
        n,rate=fixed_metric_dot_and_rate(E,K,JJ,P,gg,1.,W)
        self.assertLess(np.max(abs(np.gradient(n,x)[1:-1]-rate[1:-1])),.002)
    def test_pullback_zero_crossings_can_differ(self):
        x=np.linspace(1,8,281)
        lab=(x-2.4)*(x-4.3)*(x-6.6)
        inv=np.ones_like(x)
        h=np.gradient(lab,x);zero=np.zeros_like(x)
        self.assertEqual(len(roots_from_samples(x,lab,h)),3)
        self.assertEqual(len(roots_from_samples(x,inv,zero)),0)
        v,d=mixed_metric(lab,h,inv,zero,1.)
        self.assertEqual(len(roots_from_samples(x,v,d)),0)
    def test_dynamic_metric_rate_cancellation(self):
        x=np.linspace(0,1,5); j0=np.diag(np.linspace(.1,.6,11));j=np.tile(j0,(5,1,1))
        P=np.array([np.diag(np.exp(np.diag(j0)*t)) for t in x]);B=np.linalg.inv(P)
        W=np.einsum('nji,njk->nik',B,B)
        wp=-np.einsum('nji,njk->nik',j,W)-np.einsum('nij,njk->nik',W,j)
        early=np.ones(11);g=np.tile(np.arange(1,12)/50,(5,1));total=np.array([early+t*g[0] for t in x])
        E=np.einsum('nij,j->ni',P,early);K=np.einsum('nij,nj->ni',P,total)
        rate=dot_derivative_variable_metric(E,K,j,P,g,1.,W,wp)
        expected=np.einsum('j,nj->n',early,g)
        np.testing.assert_allclose(rate,expected,rtol=0,atol=1e-11)
    def test_reference_metric_spd(self):
        P=np.diag(np.linspace(.5,5.,11));W=fixed_w_from_reference(P)
        self.assertGreater(np.linalg.eigvalsh(W).min(),0)
        np.testing.assert_allclose(P.T@W@P,np.eye(11),atol=1e-12)

if __name__=='__main__':unittest.main()
