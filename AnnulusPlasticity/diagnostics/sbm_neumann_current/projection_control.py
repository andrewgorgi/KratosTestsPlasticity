"""Diagnostic: correct inner projection coordinates before condition initialization.

Uses the same installed condition. Only projection-node positions/normals change;
production geometries and modeler sources are untouched.
"""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker

output = HERE/'exact_inner_projection'
output.mkdir(exist_ok=True)
for filename in ('build.json', 'parameters_input.json', 'materials_input.json'):
    (output/filename).write_bytes((HERE/filename).read_bytes())
Analysis = worker.study.base.Analysis
original = Analysis.ModifyInitialGeometry


def modify(self):
    original(self)
    part = self.model['IgaModelPart']
    skin = self.model['skin_model_part'].GetSubModelPart('inner')
    def collect(part):
        nodes = list(part.Nodes) if part.Name=='projection_data' else []
        for child in part.SubModelParts:
            nodes += collect(child)
        return nodes
    nodes = sorted(collect(skin), key=lambda n:n.Id)
    conditions = sorted([c for c in part.Conditions if 'SbmLoadSolidCondition' in worker.study.base.condition_name(c)], key=lambda c:c.Id)
    assert len(nodes)==len(conditions) and nodes
    shifts = []
    for node, condition in zip(nodes, conditions):
        center = condition.GetGeometry().Center()
        exact = np.array([center.X, center.Y]); exact /= np.linalg.norm(exact)
        shifts.append(float(np.linalg.norm(exact-[node.X0, node.Y0])))
        node.X, node.Y, node.X0, node.Y0 = float(exact[0]), float(exact[1]), float(exact[0]), float(exact[1])
        node.SetValue(worker.study.KM.NORMAL, [-float(exact[0]), -float(exact[1]), 0.])
    (output/'projection_changes.json').write_text(json.dumps(dict(nodes=len(nodes), max_shift=max(shifts)), indent=2)+'\n')


Analysis.ModifyInitialGeometry = modify
sys.argv = [__file__, '--mode', 'new', '--degree', '3', '--spans', '20', '--output', str(output)]
worker.main()
