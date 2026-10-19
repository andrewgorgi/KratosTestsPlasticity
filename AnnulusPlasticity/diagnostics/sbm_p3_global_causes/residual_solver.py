"""Diagnostic Newton: same native linear solver, Armijo on the residual norm.

Keep standard solution-step initialize/finalize and all conditions unchanged.
Use stricter absolute equilibrium/correction checks, never accept failed states.
"""
import json

import numpy as np
from KratosMultiphysics.scipy_conversion_tools import to_csr


def install(study, output):
    KM = study.KM
    initialize = study.base.Analysis.Initialize
    output.mkdir(parents=True, exist_ok=True)

    def initialize_and_wrap(self):
        initialize(self)
        solver = self._GetSolver()
        from branch_points import BoundaryPoints
        material_points = BoundaryPoints(study, self.model)

        def solve():
            strategy = solver._GetSolutionStrategy()
            builder, scheme = solver._GetBuilderAndSolver(), solver._GetScheme()
            part = solver.GetComputingModelPart()
            matrix, rhs, dx = strategy.GetSystemMatrix(), strategy.GetSystemVector(), strategy.GetSolutionVector()
            dofs = sorted([d for d in builder.GetDofSet() if d.EquationId<len(rhs)], key=lambda d:d.EquationId)
            assert len(dofs)==len(rhs) and not any(d.IsFixed() for d in dofs)
            space = KM.UblasSparseSpace()
            rows = []
            time = float(self.time)

            def save(status, norm):
                (output/f't{time:.8f}.json').write_text(json.dumps(dict(time=time, status=status,
                    residual_norm=norm, iterations=rows,
                    convergence='absolute residual norm < 1e-10 and Newton correction norm < 1e-12',
                    linear_solver='same configured Kratos solver'), indent=2)+'\n')

            def set_u(u):
                for dof, value in zip(dofs, u):
                    dof.SetSolutionStepValue(float(value))

            for iteration in range(100):
                part.ProcessInfo[KM.NL_ITERATION_NUMBER] = iteration+1
                scheme.InitializeNonLinIteration(part, matrix, dx, rhs)
                space.SetToZeroMatrix(matrix); space.SetToZeroVector(rhs); space.SetToZeroVector(dx)
                builder.Build(scheme, part, matrix, rhs)
                residual = np.array(rhs)
                norm = float(np.linalg.norm(residual))
                solver._GetLinearSolver().Solve(matrix, dx, rhs)
                delta = np.array(dx)
                delta_norm = float(np.linalg.norm(delta))
                linear_error = float(np.linalg.norm(to_csr(matrix)@delta-residual)/max(norm, 1e-30))
                if norm<1e-10 and delta_norm<1e-12:
                    scheme.FinalizeNonLinIteration(part, matrix, dx, rhs)
                    material_points.commit()
                    save('converged', norm)
                    return True
                u = np.array([d.GetSolutionStepValue() for d in dofs])
                alpha, accepted = 1., False
                trial_norm = norm
                for trial in range(30):
                    set_u(u+alpha*delta)
                    space.SetToZeroVector(rhs)
                    builder.BuildRHS(scheme, part, rhs)
                    trial_norm = float(np.linalg.norm(np.array(rhs)))
                    if trial_norm <= (1.-1e-4*alpha)*norm:
                        accepted = True
                        break
                    alpha *= .5
                rows.append(dict(iteration=iteration+1, residual_norm=norm, alpha=alpha,
                    trial_residual_norm=trial_norm, accepted=accepted, correction_norm=delta_norm,
                    linear_solve_relative_error=linear_error))
                if not accepted:
                    set_u(u)
                    (output/f'branches_t{time:.8f}.json').write_text(json.dumps(
                        material_points.check(delta), indent=2)+'\n')
                    directional = []
                    for epsilon in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
                        set_u(u+epsilon*delta)
                        space.SetToZeroVector(rhs)
                        builder.BuildRHS(scheme, part, rhs)
                        sampled = np.array(rhs)
                        actual_derivative = (sampled-residual)/epsilon
                        directional.append(dict(alpha=epsilon,
                            relative_directional_derivative_mismatch=float(np.linalg.norm(actual_derivative+residual)/max(norm, 1e-30)),
                            phi_forward_slope=float((sampled@sampled-residual@residual)/(2*epsilon)),
                            phi_predicted_slope=-norm**2))
                    set_u(u)
                    (output/f'directional_t{time:.8f}.json').write_text(json.dumps(directional, indent=2)+'\n')
                    save('stalled', norm)
                    raise RuntimeError(f'Diagnostic residual Newton stalled at t={time:g}.')
                for j, value in enumerate(alpha*delta):
                    dx[j] = float(value)
                scheme.FinalizeNonLinIteration(part, matrix, dx, rhs)
            save('iteration_limit', trial_norm)
            raise RuntimeError(f'Diagnostic residual Newton exceeded 100 iterations at t={time:g}.')

        solver.SolveSolutionStep = solve

    study.base.Analysis.Initialize = initialize_and_wrap
