"""Independent synthetic mathematical verification; NOT Q550 physical measurements."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from source_transport_math import source_hessian_via_jac,transported_integrand,simpson_cell,signed_even,relative_error,angle_deg,slope_finite_amplitude
class TestSecondOrderSource(unittest.TestCase):
    def test_quadratic_rhs_hessian(self):
        A=np.arange(1,12,dtype=float)/11
        def jac(t,y):return np.diag(3.+2.*A*y)
        y=np.linspace(-.5,.5,11);u=np.arange(1,12,dtype=float)/4
        source=source_hessian_via_jac(jac,0.,y,u,2e-4)
        self.assertTrue(np.allclose(source,A*u*u,rtol=1e-9,atol=1e-9))
    def test_linear_zero_source(self):
        mat=np.diag(np.linspace(1.,2.,11))
        self.assertLess(np.linalg.norm(source_hessian_via_jac(lambda t,y:mat,0.,np.ones(11),np.ones(11),1e-3)),1e-13)
    def test_duhamel_scalar_example(self):
        # y'=a*y², initial y0=0, tangent u=1, S2=a, transport P=1.
        a=.7;ds=.1; source=source_hessian_via_jac(lambda t,y: np.array([[2*a*y[0]]]),0.,np.array([0.]),np.array([1.]),1e-3)
        g,_=transported_integrand(np.eye(1),source)
        q=sum(simpson_cell(g,g,g,ds) for _ in range(10))
        self.assertAlmostEqual(float(q[0]),a,places=11)
        for amp in (.001,.002):
            # exact y_plus(t) = +amp/(1-a*amp*t), y_minus(t) = -amp/(1+a*amp*t)
            even=0.5*(amp/(1-a*amp)-amp/(1+a*amp))
            self.assertLess(abs(even/amp**2-a),5e-6)
    def test_nested_simpson_cubic(self):
        f=lambda s:np.array([s**3])
        total=simpson_cell(f(0),f(.5),f(1),1)
        self.assertAlmostEqual(float(total[0]),.25,places=14)
    def test_signed_even(self):
        y=np.array([2.,3.]);a=np.array([.1,.2]);b=np.array([.3,-.4]);ep=.05
        yp=y+ep*a+ep**2*b;ym=y-ep*a+ep**2*b
        self.assertTrue(np.allclose(signed_even(yp,ym,y)/ep**2,b))
    def test_extrapolation(self):
        a=np.array([1e-4,2e-4,4e-4]);k=np.arange(1.,12.)
        vals=[k+(v**2)*7*k for v in a]
        self.assertTrue(np.allclose(slope_finite_amplitude(vals,a),k,atol=2e-13))
    def test_angles(self):
        self.assertAlmostEqual(angle_deg([1,0],[0,1]),90.)
        self.assertAlmostEqual(relative_error([1.,1.],[1.,1.]),0.)
if __name__=='__main__':unittest.main()
