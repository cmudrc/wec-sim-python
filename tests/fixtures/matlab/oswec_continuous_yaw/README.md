# Repository-owned OSWEC continuous-heading source

This is a derived control for the published 250 s `PassiveYawRegression`
irregular sea. The [source workflow](../../../../.github/workflows/oswec-passive-yaw-continuous-step.yml)
checks out WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and the
unchanged WEC-Sim Applications input at
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, then applies the
[repository-owned patch](../oswec_held_yaw/continuous_heading.patch). Only
the two body yaw thresholds change from 1° to zero. Zero interpolates the
heading coefficients every step, so this control has no discrete held-heading
events. The published one-degree case remains qualified.

The six raw MATLAB R2025b `.mat` files retain the 0.01 and 0.005 s runs for
each of phase seeds 1–3. Each contains time, flap yaw and yaw speed, PTO
torque and power, incident wave time and elevation, the time step, seed, and
integrated PTO work. `manifest.json` records the exact CI run, commit, patch
hash, source-file hashes, and audited full-duration step changes. The archive
was copied from that run using
[`tools/archive_oswec_continuous_yaw.py`](../../../../tools/archive_oswec_continuous_yaw.py),
which rejects missing files, hash mismatches, or failed source limits.
The files came from [CI run 39883249936](https://github.com/cmudrc/wec-sim-python/actions/runs/39883249936)
at `8ea64f3a5c37a68ff4b740ad63c80bc459ba3a3d`; all three live MATLAB
source and independent physical Python jobs passed. Wave samples match
exactly between source steps. The audited source changes are:

| Seed | Max yaw (rad) | Max yaw speed (rad/s) | Max PTO torque (N m) | PTO work change |
| --- | ---: | ---: | ---: | ---: |
| 1 | 2.075e-5 | 4.770e-6 | 0.5724 | +0.001289% |
| 2 | 1.834e-5 | 3.660e-6 | 0.4392 | -0.006538% |
| 3 | 1.068e-5 | 1.790e-6 | 0.2147 | +0.00001065% |

The [portable gate](../../../test_oswec_continuous_yaw_portable.py) uses the
same full-duration source limits and independently advances all three public
Python physical trajectories from the pinned numerical HDF5 and generated
MATLAB phase seeds. It compares every matched sample of yaw, yaw speed, PTO
torque, and wave and compares integrated PTO work. No source states, forces,
or event times are inputs to the Python solve.
