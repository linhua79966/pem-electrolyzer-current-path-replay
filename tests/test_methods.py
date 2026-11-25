"""Only synthetic in-memory test data. No original dataset is read."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from nlr_identification import state,fit1,fit2
from juelich_identification import gram_fit,coeff
from response import zoh,matched_paths,signed_terms
from lag import selftests
class Methods(unittest.TestCase):
    def test_zoh_reset(self):
        x=[0.,10.,10.,0.,0.];z=state(x,2.,(('a',0,3),('b',3,5)))
        np.testing.assert_allclose(z,[0,0,10*(-math.expm1(-.5)),0,0])
    def test_nnls(self):
        a=[1.,0.,2.,1.];b=[0.,1.,1.,2.];y=[2*x+3*z for x,z in zip(a,b)];ix=list(range(4));v,s,_=fit2(a,b,y,ix);self.assertEqual(s,'ELIGIBLE');np.testing.assert_allclose(v,[2,3]);self.assertEqual(fit2(a,a,y,ix)[1],'COLLINEAR_DESIGN');self.assertEqual(fit2(a,b,a,ix)[1],'BOUNDARY_OR_SMALL_WEIGHT');self.assertEqual(fit1(a,[-v for v in a],ix),0)
    def test_equal_file_gram(self):
        X=np.array([[1.,0.],[0.,1.],[2.,1.],[1.,2.]]);y=X@np.array([2.,3.]);G,g=gram_fit([(X,y),(np.tile(X,(3,1)),np.tile(y,3))]);v,ratio,status=coeff(G,g,[0,1]);self.assertEqual(status,'ELIGIBLE');np.testing.assert_allclose(v,[2,3]);self.assertEqual(coeff(np.ones((2,2)),np.ones(2),[0,1])[2],'COLLINEAR')
    def test_irregular_zoh(self):
        t=[0,1,3,4];u=[100,200,300,400];a=zoh(t,u,2);self.assertAlmostEqual(a[2],math.exp(-1)*100+(1-math.exp(-1))*200);self.assertRaises(ValueError,zoh,[0,1,1],[1,2,3],2)
    def test_matched_paths(self):
        t=np.arange(50.);a=100+np.sin(t)*10;c=120+np.cos(t)*20;tau=8;beta=.007
        x=matched_paths(t,a,c,4,1,tau,beta);self.assertEqual(len(set(x.values())),1)
        x=matched_paths(t,a,c,4,30,tau,beta);self.assertAlmostEqual(x['A'],beta*zoh(t,a,tau)[34]);a2=a.copy();a2[5:]+=100;x2=matched_paths(t,a2,c,4,30,tau,beta);self.assertEqual(x['H'],x2['H']);self.assertEqual(x['U'],x2['U']);self.assertNotEqual(x['A'],x2['A'])
    def test_signed_identity(self):
        z=signed_terms([3,2,5],[1,1,2],[2,2,4],[3,1,5]);self.assertAlmostEqual(z['mse_total'],sum(v for k,v in z.items() if k!='mse_total'))
    def test_lag(self):self.assertEqual(len(selftests()),12)
if __name__=='__main__':unittest.main(verbosity=2)
