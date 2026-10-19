"""Solve cyclic mixed annular torsion using either existing macro branch."""
import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import annulus_study as study
import reference
from local import material


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('on', 'off'), required=True)
    parser.add_argument('--spans', type=int, required=True)
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--dt', type=float, default=.00025)
    args = parser.parse_args()
    study.KM.ParallelUtilities.SetNumThreads(1)
    validation = reference.verify()
    study.evaluate = reference.evaluate
    study.boundary_traction = reference.boundary_traction
    output = HERE/'global'/args.mode/f'p{args.degree}_n{args.spans}'
    output.mkdir(parents=True, exist_ok=True)
    records = []

    class Analysis(study.Analysis):
        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            s = self.samples
            u, grad, stress, alpha = reference.evaluate(s['points'], self.time)
            w = s['weights']
            part = self.model['IgaModelPart']
            # Fixed amplitude normalizations remain meaningful at zero torque.
            u_scale = reference.evaluate(s['points'], .001)[0]
            stress_scale = reference.evaluate(s['points'], .001)[2]
            info = part.ProcessInfo
            record = dict(time=float(self.time), level=float(np.interp(self.time, reference.TIMES, reference.LEVELS)),
                iteration=int(info[study.KM.NL_ITERATION_NUMBER]),
                L2_amplitude_relative=float(np.sqrt(np.sum(w[:, None]*(s['displacement']-u)**2)
                    /np.sum(w[:, None]*u_scale**2))),
                stress_amplitude_relative=float(np.sqrt(np.sum(w[:, None]*(s['stress'][:, :3]-stress)**2)
                    /np.sum(w[:, None]*stress_scale**2))),
                alpha_rms_error=float(np.sqrt(np.sum(w*(s['equivalent_plastic_strain']-alpha)**2)/w.sum())))
            records.append(record)
            np.savez(output/'last_solution.npz', **s)
            (output/'steps.json').write_text(json.dumps(records, indent=2)+'\n')
            print('CYCLIC_STEP', json.dumps(record), flush=True)

    settings = SimpleNamespace(parameters=str(HERE/'parameters_input.json'), boundary='mixed',
        degree=args.degree, dt=args.dt, echo_level=0)
    parameters = study.load_parameters(ROOT/'sbm', settings, [args.spans]*2)
    parameters['problem_data']['end_time'].SetDouble(reference.END_TIME)
    material_path = material('linear')
    variables = json.loads(material_path.read_text())['properties'][0]['Material']['Variables']
    for name, expected in {'YOUNG_MODULUS': reference.E, 'POISSON_RATIO': reference.NU,
                           'ISOTROPIC_HARDENING_MODULUS': reference.H, 'YIELD_STRESS': reference.YIELD,
                           'EXPONENTIAL_SATURATION_YIELD_STRESS': reference.YIELD,
                           'HARDENING_EXPONENT': 0.}.items():
        assert variables[name] == expected, f'Reference/material mismatch: {name}'
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(material_path))
    (output/'parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    with study.shared.working_directory(ROOT/'sbm'):
        analysis = Analysis(study.KM.Model(), parameters)
        try:
            analysis.Run()
            result = dict(status='converged', time=float(analysis.time))
        except RuntimeError as error:
            result = dict(status='failed', time=float(analysis.time), error=str(error))
    result.update(mode=args.mode, degree=args.degree, spans=args.spans, dt=args.dt,
                  reference_validation=validation, steps=records)
    (output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print('CYCLIC_RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
