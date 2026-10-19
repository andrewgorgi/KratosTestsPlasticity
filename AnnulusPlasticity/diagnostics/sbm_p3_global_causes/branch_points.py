"""Track independent boundary material histories for a loading-direction check."""
from math import factorial
import numpy as np

from check_projection_history import MaterialPoint, strain_matrix


class BoundaryPoints:
    def __init__(self, study, model):
        self.KM = study.KM
        part = model['IgaModelPart']
        conditions = sorted([c for c in part.Conditions if 'SbmLoadSolidCondition' in study.base.condition_name(c)], key=lambda c:c.Id)
        def collect(part):
            nodes = list(part.Nodes) if part.Name=='projection_data' else []
            for child in part.SubModelParts:
                nodes += collect(child)
            return nodes
        projections = sorted(collect(model['skin_model_part'].GetSubModelPart('inner')), key=lambda n:n.Id)
        assert len(conditions)==len(projections)
        self.points = []
        for c, node in zip(conditions, projections):
            g = c.GetGeometry()
            center = g.Center()
            shift = np.array([node.X0-center.X, node.Y0-center.Y])
            gradients = np.array(g.ShapeFunctionDerivatives(1, 0))
            extended = gradients.copy()
            degree = round(np.sqrt(len(g)))-1
            for order in range(2, 2*degree+1):
                d = np.array(g.ShapeFunctionDerivatives(order, 0))
                for k in range(order):
                    coefficient = shift[0]**(order-1-k)*shift[1]**k/(factorial(order-1-k)*factorial(k))
                    extended[:, 0] += coefficient*d[:, k]
                    extended[:, 1] += coefficient*d[:, k+1]
            self.points.append((c, strain_matrix(gradients), strain_matrix(extended),
                                MaterialPoint(c, part.ProcessInfo), MaterialPoint(c, part.ProcessInfo)))
        self.elastic = np.diag([67000., 67000., 33500.])

    def displacement(self, condition):
        return np.array([[n.GetSolutionStepValue(self.KM.DISPLACEMENT_X),
                          n.GetSolutionStepValue(self.KM.DISPLACEMENT_Y)] for n in condition.GetGeometry()]).ravel()

    def commit(self):
        for c, Bs, Bt, s, t in self.points:
            u = self.displacement(c)
            s.commit(Bs@u); t.commit(Bt@u)

    def check(self, direction, alpha=1e-7):
        counts = {side:dict(points=0, base_active_plastic=0, trial_active_plastic=0,
                           elastic_to_plastic=0, plastic_to_elastic=0) for side in ('surrogate', 'true')}
        norm = np.linalg.norm(self.elastic)
        for c, Bs, Bt, s, t in self.points:
            u = self.displacement(c)
            local_direction = np.array([direction[d.EquationId] for d in c.GetDofList(s.values.GetProcessInfo())])
            for name, B, point in [('surrogate', Bs, s), ('true', Bt, t)]:
                _, D0 = point.response(B@u)
                _, D1 = point.response(B@(u+alpha*local_direction))
                active0 = np.linalg.norm(D0-self.elastic)/norm > 1e-8
                active1 = np.linalg.norm(D1-self.elastic)/norm > 1e-8
                row = counts[name]
                row['points'] += 1
                row['base_active_plastic'] += int(active0)
                row['trial_active_plastic'] += int(active1)
                row['elastic_to_plastic'] += int(not active0 and active1)
                row['plastic_to_elastic'] += int(active0 and not active1)
        return dict(alpha=alpha, counts=counts)
