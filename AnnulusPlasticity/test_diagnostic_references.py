"""Validate the independent perfect-plastic manufactured controls."""
import unittest
import numpy as np
from diagnose_mixed import reference
from analytic import YIELD


class DiagnosticReferenceTests(unittest.TestCase):
    def test_perfect_patch_is_uniform_plastic_with_zero_body_force(self):
        evaluate,body = reference('patch_perfect')
        points = np.array([[1.,0.],[0.,1.],[-1.5,.2]])
        u,gradient,stress = evaluate(points,.002)
        np.testing.assert_allclose(u[:,0],.02*points[:,0])
        np.testing.assert_allclose(gradient,np.tile([.02,0,0,0],(len(points),1)))
        np.testing.assert_allclose(body(points,.002),0)
        xx,yy,xy,zz = stress.T
        q = np.sqrt(((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)/2+3*xy**2)
        np.testing.assert_allclose(q,YIELD,rtol=1e-12)

    def test_nonaffine_perfect_body_force_matches_stress_divergence(self):
        evaluate,body = reference('cosine_perfect')
        points = np.array([[1.3,.4],[-.8,1.3],[.4,-1.5]])
        for time in (1e-5,.002):
            h = 1e-6
            derivatives = []
            for axis in range(2):
                delta = np.eye(2)[axis]*h
                derivatives.append((evaluate(points+delta,time)[2]-evaluate(points-delta,time)[2])/(2*h))
            expected = np.column_stack((-(derivatives[0][:,0]+derivatives[1][:,2]),
                                         -(derivatives[0][:,2]+derivatives[1][:,1]),
                                         np.zeros(len(points))))
            np.testing.assert_allclose(body(points,time),expected,rtol=1e-6,atol=1e-6)


if __name__ == '__main__':
    unittest.main()
