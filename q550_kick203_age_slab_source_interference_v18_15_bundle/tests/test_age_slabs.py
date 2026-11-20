import sys,unittest
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from age_slab_math import slab_edges,metrics,orth,slabs

class TestAgeSlabs(unittest.TestCase):
    def setUp(self):
        self.x=np.linspace(0,8,641)
        self.e=np.array([1.,0.,0.,0.])
        # Synthetic history with alternating source sign and independent lateral motion.
        self.A=np.stack([np.zeros_like(self.x),self.x/5+np.sin(2*np.pi*self.x),
                        np.cos(np.pi*self.x)-1, .3*self.x],axis=-1)
        self.ai=PchipInterpolator(self.x,self.A,axis=0)
    def test_partition_covers_exact_age(self):
        for h in (.5,.25,.125,.0625):
            edges=slab_edges(3.9496,h)
            self.assertAlmostEqual(edges[0],0)
            self.assertAlmostEqual(edges[-1],3.9496)
            self.assertTrue(np.all(np.diff(edges)>0))
    def test_vector_and_pairwise_closure(self):
        for h in (.5,.25,.125,.0625):
            q,rows,_=metrics(self.ai,self.e,3.9496,h)
            self.assertLess(q['closure_error'],1e-11)
            self.assertLess(q['pair_identity_error'],1e-11)
            self.assertLess(q['axial_closure_error'],1e-11)
            self.assertAlmostEqual(sum(r['signed_axial'] for r in rows),q['rho'])
    def test_negative_source_slab_can_oppose(self):
        r,rows,_=metrics(self.ai,self.e,3.9496,.125)
        self.assertGreater(r['negative_axial'],0)
        self.assertTrue(any(row['role']=='opposing' for row in rows))
    def test_full_sum_independent_of_resolution(self):
        a=[metrics(self.ai,self.e,4.123,h)[0]['rho'] for h in (.5,.25,.125,.0625)]
        self.assertLess(max(a)-min(a),1e-12)
    def test_coordinates_preserve_linear_identity(self):
        T=np.array([[2,.3,0,0],[0,1,.4,.1],[0,0,1,.25],[0,0,0,.75]])
        e=self.e;edges,s=slabs(self.ai,e,3.9496,.25)
        orig=np.stack([e,s[0],s[1]],axis=1)
        self.assertEqual(np.linalg.matrix_rank(orig),np.linalg.matrix_rank(T@orig))
    def test_metric_dependence_of_slabs(self):
        r,rows,_=metrics(self.ai,self.e,3.9496,.25)
        T=np.diag([1.,2.,.7,1.3]);a=PchipInterpolator(self.x,self.A@T.T,axis=0)
        rr,rows2,_=metrics(a,T@self.e,3.9496,.25)
        self.assertGreater(abs(r['negative_axial']-rr['negative_axial']),1e-4)
    def test_zero_vector_rejected(self):
        const=PchipInterpolator(self.x,np.tile(self.e,(len(self.x),1)),axis=0)
        with self.assertRaises(ValueError):metrics(const,self.e,3.,.25)
    def test_old2_recent2_split(self):
        z=3.9496
        r,rows,_=metrics(self.ai,self.e,z,.125)
        recent=orth(self.e,self.ai(z)-self.ai(z-2));old=orth(self.e,self.ai(z-2))
        self.assertAlmostEqual(r['old2_norm'],np.linalg.norm(old),places=10)
        self.assertAlmostEqual(r['recent2_norm'],np.linalg.norm(recent),places=10)

if __name__=='__main__':unittest.main()
