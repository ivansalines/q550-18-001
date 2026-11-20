import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from echo_math import compensate,transport_from_origin,signed_even_odd,safe_ratio,unit_angle

class TestEcho(unittest.TestCase):
    def test_linear_zero_lag(self):
        y=np.array([3.,2.]);v=np.array([1.,-.2]);e=.003
        for sig in (-1,1):
            r=compensate(y+sig*e*v,y,v,e,sig)
            self.assertLess(np.linalg.norm(r-y),1e-14)
    def test_linear_evolved_echo(self):
        Phi=np.diag([2.,.7]); y0=np.zeros(2);v=np.array([.8,.6]);a=.005
        for sig in (-1,1):
            write=Phi@(y0+sig*a*v)
            reverse=compensate(write,Phi@y0,Phi@v,a,sig)
            self.assertLess(np.linalg.norm(reverse-Phi@y0),1e-14)
    def test_nonlinear_quadratic_echo(self):
        # Explicit forward flow of dy/ds=y+y^2 with baseline y=0.
        def flow(x,t):
            e=np.exp(t)
            return x*e/(1+x*(1-e))
        t=.5
        result=[]
        for amp in (.002,.001):
            residual=[]
            for sig in (-1,1):
                residual.append(flow(sig*amp,t)-sig*amp*np.exp(t))
            even,odd=signed_even_odd(residual[1],residual[0])
            result.append((abs(even),abs(odd)))
        # Doubling epsilon gives approximately 4x even, ~8x odd.
        self.assertTrue(3.9<result[0][0]/result[1][0]<4.1)
        self.assertTrue(7.8<result[0][1]/result[1][1]<8.2)
    def test_transport(self):
        P0=np.diag([2.,3.]);P1=np.diag([4.,12.]);r=np.array([1.,2.])
        pr,c=transport_from_origin(P1,P0,r)
        np.testing.assert_allclose(pr,[2.,8.]);self.assertAlmostEqual(c,1.5)
    def test_condition_gate(self):
        P=np.diag([1.,1e-12]);a,c=transport_from_origin(np.eye(2),P,np.array([.1,.2]),1e10)
        self.assertTrue(np.isnan(a).all());self.assertGreater(c,1e10)
    def test_even_odd_and_angle(self):
        e,o=signed_even_odd(np.array([3.,2.]),np.array([1.,4.]))
        np.testing.assert_allclose(e,[2.,3.]);np.testing.assert_allclose(o,[1.,-1.]);self.assertAlmostEqual(unit_angle(e,e),0,places=4)

if __name__=='__main__':unittest.main()
