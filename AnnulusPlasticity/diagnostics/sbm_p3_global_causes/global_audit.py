"""Global spectra/energy attribution and a diagnostic residual-norm Newton solve."""
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvals, eigh, svdvals
from scipy.sparse import save_npz
from scipy.sparse.linalg import spsolve
from KratosMultiphysics.scipy_conversion_tools import to_csr


def install(study, output):
    KM = study.KM
    output.mkdir(parents=True, exist_ok=True)
    initialize = study.base.Analysis.Initialize

    def initialize_and_wrap(self):
        initialize(self)
        solver = self._GetSolver()
        native = solver.SolveSolutionStep

        def context():
            strategy = solver._GetSolutionStrategy()
            return (solver._GetBuilderAndSolver(), solver._GetScheme(),
                    solver.GetComputingModelPart(), strategy.GetSystemMatrix(), strategy.GetSystemVector())

        def build():
            builder, scheme, part, matrix, rhs = context()
            space = KM.UblasSparseSpace()
            space.SetToZeroMatrix(matrix); space.SetToZeroVector(rhs)
            builder.Build(scheme, part, matrix, rhs)
            return to_csr(matrix).copy(), np.array(rhs)

        def active_dofs():
            builder, _, _, _, rhs = context()
            return sorted([d for d in builder.GetDofSet() if d.EquationId < len(rhs)], key=lambda d:d.EquationId)

        def vector():
            return np.array([d.GetSolutionStepValue() for d in active_dofs()])

        def set_vector(u):
            for dof, value in zip(active_dofs(), u):
                dof.SetSolutionStepValue(float(value))

        def residual():
            builder, scheme, part, _, rhs = context()
            KM.UblasSparseSpace().SetToZeroVector(rhs)
            builder.BuildRHS(scheme, part, rhs)
            return np.array(rhs)

        def snapshot(name):
            folder = output/name
            folder.mkdir(exist_ok=True)
            A, R = build()
            save_npz(folder/'matrix.npz', A)
            dense = A.toarray()
            eig = eigvals(dense)
            symmetric, modes = eigh((dense+dense.T)/2, subset_by_index=(0, 0))
            singular = svdvals(dense)
            mode = modes[:, 0]
            energies = defaultdict(float)
            _, _, part, _, _ = context()
            for entities in (part.Elements, part.Conditions):
                for entity in entities:
                    if entity.IsDefined(KM.ACTIVE) and not entity.Is(KM.ACTIVE):
                        continue
                    local_dofs = entity.GetDofList(part.ProcessInfo)
                    if not len(local_dofs):
                        continue
                    local = np.array([mode[d.EquationId] if d.EquationId < len(mode) else 0. for d in local_dofs])
                    K, b = KM.Matrix(), KM.Vector()
                    entity.CalculateLocalSystem(K, b, part.ProcessInfo)
                    energies[entity.Info().split('#')[0].strip()] += float(local@np.array(K)@local)
            dofs = active_dofs()
            coordinates = np.array([[part.GetNode(d.Id()).X0, part.GetNode(d.Id()).Y0] for d in dofs])
            radii = np.linalg.norm(coordinates, axis=1)
            row = dict(time=float(self.time), residual_norm=float(np.linalg.norm(R)), dofs=len(R),
                min_real_eigenvalue=float(eig.real.min()), negative_real_eigenvalues=int((eig.real < -1e-9*max(abs(eig))).sum()),
                symmetric_min_eigenvalue=float(symmetric[0]), smallest_singular_value=float(singular[-1]),
                condition_number_2=float(singular[0]/singular[-1]), mode_energy_by_entity=dict(energies),
                mode_energy_reconstruction_error=float(abs(sum(energies.values())-symmetric[0])/max(1., abs(symmetric[0]))),
                mode_amplitude_fraction_r_less_1p3=float(np.sum(mode[radii<1.3]**2)/np.sum(mode**2)))
            (folder/'summary.json').write_text(json.dumps(row, indent=2)+'\n')
            np.savez(folder/'mode.npz', coordinates=coordinates, mode=mode)
            return row

        def armijo(start):
            failed = vector()
            set_vector(start)
            rows = []
            outcome = 'iteration_limit'
            try:
                for iteration in range(100):
                    A, R = build()
                    norm = float(np.linalg.norm(R))
                    if norm < 1e-10:
                        outcome = 'residual_converged'
                        break
                    delta = spsolve(A, R)
                    linear_error = float(np.linalg.norm(A@delta-R)/max(norm, 1e-30))
                    current = vector()
                    alpha = 1.
                    accepted = False
                    trial_norm = norm
                    for trial in range(30):
                        set_vector(current+alpha*delta)
                        trial_norm = float(np.linalg.norm(residual()))
                        if trial_norm <= (1.-1e-4*alpha)*norm:
                            accepted = True
                            break
                        alpha *= .5
                    rows.append(dict(iteration=iteration, residual_norm=norm, alpha=alpha,
                        trial_residual_norm=trial_norm, accepted=accepted, linear_solve_relative_error=linear_error,
                        displacement_norm=float(np.linalg.norm(current)), increment_norm=float(np.linalg.norm(delta))))
                    if not accepted:
                        set_vector(current)
                        outcome = 'line_search_stalled'
                        break
                _, R = build()
                data = dict(status=outcome, time=float(self.time), final_residual_norm=float(np.linalg.norm(R)),
                            iterations=rows, committed=False)
                (output/'residual_armijo.json').write_text(json.dumps(data, indent=2)+'\n')
            finally:
                set_vector(failed)
                build()

        def checked():
            # At the first solve the DOF set may not yet have been constructed.
            previous = vector() if len(solver._GetBuilderAndSolver().GetDofSet()) else None
            try:
                converged = native()
            except RuntimeError:
                snapshot('failed')
                from audit_global_cycle import audit
                audit(solver, output/'failed')
                if previous is not None:
                    armijo(previous)
                raise
            if any(abs(float(self.time)-t)<1e-12 for t in (.00025, .0005)):
                snapshot('elastic' if float(self.time)<.0003 else 'plastic')
            return converged

        solver.SolveSolutionStep = checked

    study.base.Analysis.Initialize = initialize_and_wrap
