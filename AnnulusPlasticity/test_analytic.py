"""Independent checks of the exact field, without importing Kratos."""
import unittest
import json
from pathlib import Path
import numpy as np
from analytic import evaluate, traction, boundary_traction, END_TIME, C_FINAL, G, H, YIELD


class AnalyticTests(unittest.TestCase):
    def test_inner_traction_has_material_outward_orientation(self):
        angle = np.linspace(0,2*np.pi,33)
        er = np.column_stack((np.cos(angle),np.sin(angle)))
        et = np.column_stack((-er[:,1],er[:,0]))
        for time in (0.,END_TIME/2,END_TIME):
            c = C_FINAL*time/END_TIME
            np.testing.assert_allclose(boundary_traction(er,time),-c*et,atol=2e-13)
            np.testing.assert_allclose(boundary_traction(2*er,time),c/4*et,atol=2e-13)

    def test_boundaries_and_yield(self):
        u, _, stress, alpha = evaluate([[1.,0.],[2.,0.]],END_TIME)
        np.testing.assert_allclose(u[:,1],[-0.010565772032267568,0],atol=1e-14)
        np.testing.assert_allclose(stress[:,2],[C_FINAL,C_FINAL/4])
        np.testing.assert_allclose(alpha,[.025125,0],atol=1e-14)
        np.testing.assert_allclose(traction([[1.,0.],[2.,0.]],[[-1.,0.],[1.,0.]],END_TIME),
                                   [[0.,-C_FINAL],[0.,C_FINAL/4]])

    def test_gradient_and_equilibrium(self):
        points = np.array([[1.1,.2],[.3,1.8],[-1.2,-.4]])
        h = 1e-6
        for time in (0.,END_TIME/10,END_TIME):
            u, grad, stress, alpha = evaluate(points,time)
            derivatives = []
            stress_derivatives = []
            for axis in range(2):
                delta = np.eye(2)[axis]*h
                plus, minus = evaluate(points+delta,time),evaluate(points-delta,time)
                derivatives.append((plus[0]-minus[0])/(2*h))
                stress_derivatives.append((plus[2]-minus[2])/(2*h))
            expected_grad = np.stack(derivatives,axis=2).reshape(-1,4)
            np.testing.assert_allclose(grad,expected_grad,atol=1e-10,rtol=1e-7)
            divergence = np.column_stack((stress_derivatives[0][:,0]+stress_derivatives[1][:,2],
                                          stress_derivatives[0][:,2]+stress_derivatives[1][:,1]))
            np.testing.assert_allclose(divergence,0,atol=1e-7)
            # Reconstruct radial shear from the Cartesian strain tensor.
            r = np.linalg.norm(points,axis=1)
            er = points/r[:,None]
            et = np.column_stack((-er[:,1],er[:,0]))
            eps = (grad.reshape(-1,2,2)+grad.reshape(-1,2,2).transpose(0,2,1))/2
            gamma = 2*np.einsum('ni,nij,nj->n',er,eps,et)
            tau = C_FINAL*time/END_TIME/r**2
            np.testing.assert_allclose(gamma,tau/G+np.sqrt(3)*alpha,atol=1e-14)
            plastic = alpha > 0
            np.testing.assert_allclose(np.sqrt(3)*tau[plastic],YIELD+H*alpha[plastic])

    def test_supplied_nurbs_circles(self):
        # Rational quadratic Bezier arcs: each pair of knot spans shares an endpoint.
        root = Path(__file__).resolve().parents[1]/'AnalyticalTestRevised/geometry'
        for side,radius in [('inner',1.),('outer',2.)]:
            curve = json.loads((root/f'{side}_circle.json').read_text())['Lines'][0]
            cp = np.array(curve['CPCoordinates'])
            weights = np.array(curve['Weights'])
            for start in (0,2,4):
                s = np.linspace(0,1,51)
                basis = np.column_stack(((1-s)**2,2*s*(1-s),s*s))*weights[start:start+3]
                points = basis@cp[start:start+3]/basis.sum(axis=1)[:,None]
                np.testing.assert_allclose(np.linalg.norm(points[:,:2],axis=1),radius,atol=1e-14)


if __name__ == '__main__':
    unittest.main()
