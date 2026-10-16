"""Reproduce a mixed SBM failure and audit its assembled trial-state Jacobian.

The failed step is not finalized. Perturbations preserve committed material
history. All production parameters/materials remain unchanged.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse import save_npz
from scipy.sparse.linalg import eigs
from KratosMultiphysics.scipy_conversion_tools import to_csr

import annulus_study as study

KM = study.KM


def audit(solver, output):
    builder = solver._GetBuilderAndSolver()
    scheme = solver._GetScheme()
    strategy = solver._GetSolutionStrategy()
    part = solver.GetComputingModelPart()
    matrix = strategy.GetSystemMatrix()
    rhs = strategy.GetSystemVector()
    space = KM.UblasSparseSpace()
    dofs = list(builder.GetDofSet())
    base = np.array([d.GetSolutionStepValue() for d in dofs])
    size = len(rhs)

    def assemble():
        space.SetToZeroMatrix(matrix)
        space.SetToZeroVector(rhs)
        builder.Build(scheme, part, matrix, rhs)

    def restore():
        for dof, value in zip(dofs, base):
            dof.SetSolutionStepValue(float(value))

    assemble()
    result = dict(time=float(part.ProcessInfo[KM.TIME]), dofs=size,
                  iteration=int(part.ProcessInfo[KM.NL_ITERATION_NUMBER]),
                  residual_norm=float(np.linalg.norm(np.array(rhs))), jacobian=[])
    cycle_matrix = to_csr(matrix).copy()
    save_npz(output/'cycle_matrix.npz', cycle_matrix)
    # Three independent directions; central differences at fixed old history.
    rng = np.random.default_rng(251)
    try:
        for index in range(3):
            direction = rng.normal(size=size)
            direction /= np.linalg.norm(direction)
            product = KM.Vector(size)
            space.Mult(matrix, KM.Vector(direction), product)
            expected = np.array(product)
            for h in (1e-7, 1e-8, 1e-9):
                values = []
                for sign in (1, -1):
                    for dof, value in zip(dofs, base):
                        if dof.EquationId < size:
                            dof.SetSolutionStepValue(float(
                                value + sign*h*direction[dof.EquationId]))
                    space.SetToZeroVector(rhs)
                    builder.BuildRHS(scheme, part, rhs)
                    values.append(np.array(rhs))
                fd = -(values[0]-values[1])/(2*h)
                result['jacobian'].append(dict(
                    direction=index, h=h,
                    error=float(np.linalg.norm(expected-fd)/np.linalg.norm(fd))))
            restore()
        (output/'audit.json').write_text(json.dumps(result, indent=2)+'\n')

        # First-step zero-displacement comparison: identical old history,
        # mesh, frozen Dirichlet operator and prescribed loads. It is an
        # elastic tangent comparison, not another equilibrium solution.
        if int(part.ProcessInfo[KM.STEP]) == 1:
            for dof in dofs:
                dof.SetSolutionStepValue(0.)
            assemble()
            elastic_matrix = to_csr(matrix).copy()
            save_npz(output/'elastic_matrix.npz', elastic_matrix)
            spectra = {}
            for name, operator in [('cycle', cycle_matrix), ('elastic', elastic_matrix)]:
                values = eigs(operator, k=8, sigma=0, which='LM',
                              return_eigenvectors=False, tol=1e-10)
                spectra[name] = [[float(v.real), float(v.imag)] for v in values]
            (output/'spectrum.json').write_text(json.dumps(spectra, indent=2)+'\n')
        print('GLOBAL_AUDIT', json.dumps(result), flush=True)
    finally:
        restore()
        # Return the system to the failed trial state, without committing it.
        assemble()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--spans', type=int, default=40)
    parser.add_argument('--dt', type=float, default=.001)
    parser.add_argument('--max-iterations', type=int, default=100)
    parser.add_argument('--output', type=Path,
                        default=study.ROOT/'diagnostics/global_cycle_audit')
    args = parser.parse_args()
    if min(args.degree, args.spans, args.max_iterations) < 1 or args.dt <= 0:
        parser.error('Degree, spans, iteration limit and dt must be positive.')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    KM.ParallelUtilities.SetNumThreads(1)

    class Analysis(study.Analysis):
        def Initialize(self):
            super().Initialize()
            solver = self._GetSolver()
            solve = solver.SolveSolutionStep

            def solve_and_audit():
                try:
                    return solve()
                except RuntimeError:
                    audit(solver, output)
                    raise

            solver.SolveSolutionStep = solve_and_audit

    settings = SimpleNamespace(parameters=None, boundary='mixed',
                               degree=args.degree, dt=args.dt, echo_level=1)
    folder = study.ROOT/'sbm'
    parameters = study.load_parameters(folder, settings, [args.spans]*2)
    parameters['solver_settings']['max_iteration'].SetInt(args.max_iterations)
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(
        str(study.ROOT/'materials.json'))
    (output/'parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    with study.shared.working_directory(folder):
        analysis = Analysis(KM.Model(), parameters)
        try:
            analysis.Run()
            result = dict(status='converged', time=float(analysis.time))
        except RuntimeError as error:
            result = dict(status='failed', time=float(analysis.time), error=str(error))
        (output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
        print('SOLVE_RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
