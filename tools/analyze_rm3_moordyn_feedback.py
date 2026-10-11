"""Attribute RM3 visualization drift using source-derived mass feedback.

This is a diagnostic counterfactual, not an independent parity run. The
prescribed correction comes from MATLAB's saved force and acceleration logs.
"""

import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np

from wecsim import MoorDyn
from wecsim.generalDynamics import GeneralizedDynamics
from wecsim.rm3Regular import solve_rm3_regular


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--prefix', default='dense',
                    help='prefix of the 10 s source CSV files')
parser.add_argument('records', type=Path, help='MATLAB dense CSV directory')
parser.add_argument('hydro', type=Path, help='pinned RM3 hydrodynamic HDF5')
parser.add_argument('lines', type=Path, help='pinned MoorDyn lines.txt')
parser.add_argument('library', type=Path, help='pinned native MoorDyn library')
parser.add_argument('output', type=Path, help='directory for summary and traces')
args = parser.parse_args()
records = args.records.resolve(strict=True)
hydro = args.hydro.resolve(strict=True)
lines_source = args.lines.resolve(strict=True)
library = args.library.resolve(strict=True)
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
source = [np.loadtxt(records / f'{args.prefix}_body{n}.csv', delimiter=',')
          for n in (1, 2)]
terms = [np.loadtxt(records / f'{args.prefix}_forces_body{n}.csv', delimiter=',')
         for n in (1, 2)]
assert all(values.shape == (1001, 25) for values in source)
assert all(values.shape == (1001, 49) for values in terms)
excitation = np.stack([term[:, 7:13] for term in terms], axis=1)
correction = np.zeros((1001, 4))

with h5py.File(hydro) as h5:
    for number in (1, 2):
        group = h5[f'body{number}']
        full = 1000 * np.asarray(group['hydro_coeffs/added_mass/inf_freq'])
        added = full[:, 6 * (number - 1):6 * number]
        applied_added = added.copy()
        shift = 2 * np.trace(added[:3, :3])
        applied_added[:3, :3] -= np.eye(3) * shift
        applied_added[4, 4] = 0
        center_z = float(np.asarray(group['properties/cg']).ravel()[2])
        for index in range(1001):
            angle = source[number - 1][index, 5]
            sine, cosine = np.sin(angle), np.cos(angle)
            slide = ((source[number - 1][index, 3] - center_z
                      - center_z * (cosine - 1)) / cosine)
            radius = center_z + slide
            jacobian = np.zeros((6, 4))
            jacobian[0, 0] = 1
            jacobian[0, number] = sine
            jacobian[0, 3] = radius * cosine
            jacobian[2, number] = cosine
            jacobian[2, 3] = -radius * sine
            jacobian[4, 3] = 1
            acceleration = terms[number - 1][index, 1:7]
            applied_log = terms[number - 1][index, 19:25].copy()
            applied_log[4] -= added[4, 4] * acceleration[4]
            correction[index] += jacobian.T @ (
                applied_added @ acceleration - applied_log
            )

print('feedback correction maximum by coordinate:',
      np.max(np.abs(correction), axis=0))
print('feedback correction mean by coordinate:', np.mean(correction, axis=0))

original = GeneralizedDynamics.acceleration
summary = {}
for horizon in (0.0, 0.1, 0.5, 1.0, 10.0):
    directory = output / f'horizon_{horizon:g}'
    directory.mkdir(exist_ok=True)
    lines = directory / 'lines.txt'
    shutil.copyfile(lines_source, lines)

    def augmented(self, at_time, coordinate, speed, **kwargs):
        index = round(at_time / 0.01)
        if horizon > 0 and index <= round(horizon / 0.01):
            external = kwargs.get('applied_force')
            kwargs['applied_force'] = (correction[index].copy()
                                       if external is None else
                                       external + correction[index])
        return original(self, at_time, coordinate, speed, **kwargs)

    GeneralizedDynamics.acceleration = augmented
    try:
        result = solve_rm3_regular(
            hydro, wave_height=0, wave_period=8,
            pitch_inertias=(21_306_090.66, 94_407_091.24),
            pto_damping=1_200_000, moordyn=MoorDyn(library, lines),
            moordyn_point=(0, 0, 21.5),
            radiation_memory=60, radiation_method='convolution',
            added_mass_scheme='implicit', excitation_force=excitation,
            initial_coordinate=np.array([0, 0, -0.21, 0]),
            dt=0.01, end_time=10, ramp_time=0,
        )
    finally:
        GeneralizedDynamics.acceleration = original
    for number in (0, 1):
        np.savetxt(directory / f'body{number + 1}.csv',
                   np.column_stack((result.time,
                                    result.body_position[:, number],
                                    result.body_velocity[:, number])),
                   delimiter=',')
    print('horizon', horizon)
    summary[str(horizon)] = {}
    for number, label in ((0, 'float'), (1, 'spar')):
        error = result.body_position[:, number, 0] - source[number][:, 1]
        summary[str(horizon)][label] = {
            'max_surge_error_m': float(np.max(np.abs(error))),
            'error_at_0p1s_m': float(error[10]),
            'error_at_1s_m': float(error[100]),
            'error_at_10s_m': float(error[-1]),
        }
        print(label, summary[str(horizon)][label])

(output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
