import json
from collections import defaultdict
import numpy as np
import annulus_study as study

KM = study.KM


def audit_entities(solver, output):
    part = solver.GetComputingModelPart()
    info = part.ProcessInfo
    dofs = list(solver._GetBuilderAndSolver().GetDofSet())
    size = len(solver._GetSolutionStrategy().GetSystemVector())
    base = [d.GetSolutionStepValue() for d in dofs]
    rng = np.random.default_rng(251)
    direction = rng.normal(size=size); direction /= np.linalg.norm(direction)
    entities = []
    sums = defaultdict(lambda: dict(expected=np.zeros(size), fd=np.zeros(size),
                                   count=0, max_local_error=0., worst_entity=None))
    for collection in (part.Elements, part.Conditions):
        for entity in collection:
            if entity.IsDefined(KM.ACTIVE) and not entity.Is(KM.ACTIVE):
                continue
            label = entity.Info().split('#')[0].strip()
            local_dofs = entity.GetDofList(info)
            if not len(local_dofs):
                continue
            ids = np.array([d.EquationId for d in local_dofs], dtype=np.int64)
            v = np.array([direction[i] if i < size else 0. for i in ids])
            matrix, rhs = KM.Matrix(), KM.Vector()
            entity.CalculateLocalSystem(matrix, rhs, info)
            K, R = np.array(matrix), np.array(rhs)
            product = K@v
            valid = ids < size
            np.add.at(sums[label]['expected'], ids[valid], product[valid])
            entities.append((entity, label, ids, product, R))
    h = 1e-8
    try:
        trials = []
        for sign in (1., -1.):
            for dof, value in zip(dofs, base):
                if dof.EquationId < size:
                    dof.SetSolutionStepValue(float(value+sign*h*direction[dof.EquationId]))
            values = []
            for entity, label, ids, product, original in entities:
                rhs = KM.Vector()
                entity.CalculateRightHandSide(rhs, info)
                values.append(np.array(rhs))
            trials.append(values)
    finally:
        for dof, value in zip(dofs, base):
            dof.SetSolutionStepValue(float(value))
    repeat_error = 0.
    for index, (entity, label, ids, product, original) in enumerate(entities):
        fd = -(trials[0][index]-trials[1][index])/(2*h)
        valid = ids < size
        np.add.at(sums[label]['fd'], ids[valid], fd[valid])
        entry = sums[label]
        error = float(np.linalg.norm(product-fd)/max(1., np.linalg.norm(fd)))
        entry['count'] += 1
        if error > entry['max_local_error']:
            entry['max_local_error'], entry['worst_entity'] = error, entity.Id
        rhs = KM.Vector()
        entity.CalculateRightHandSide(rhs, info)
        repeat_error = max(repeat_error, float(np.linalg.norm(np.array(rhs)-original)/max(1.,np.linalg.norm(original))))
    total_expected = sum(v['expected'] for v in sums.values())
    total_fd = sum(v['fd'] for v in sums.values())
    matrix_product = KM.Vector(size)
    KM.UblasSparseSpace().Mult(solver._GetSolutionStrategy().GetSystemMatrix(),
                              KM.Vector(direction), matrix_product)
    result = dict(repeat_error=repeat_error, by_type={},
        global_error=float(np.linalg.norm(total_expected-total_fd)/np.linalg.norm(total_fd)),
        matrix_reconstruction_error=float(np.linalg.norm(total_expected-np.array(matrix_product))
                                         /np.linalg.norm(np.array(matrix_product))))
    for label, entry in sums.items():
        expected, fd = entry.pop('expected'), entry.pop('fd')
        result['by_type'][label] = dict(**entry,
            assembled_directional_error=float(np.linalg.norm(expected-fd)/max(1.,np.linalg.norm(fd))),
            mismatch_norm=float(np.linalg.norm(expected-fd)), fd_norm=float(np.linalg.norm(fd)))
    (output/'entity_audit.json').write_text(json.dumps(result, indent=2))
    print('ENTITY_AUDIT', json.dumps(result), flush=True)
