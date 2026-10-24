# Ellipsoid no-active-delay MaxStep continuation

These fixtures archive the four successful 150 s source controls from
[GitHub Actions run 40131806849](https://github.com/cmudrc/wec-sim-python/actions/runs/40131806849)
at repository commit `f96a2056637f9bf0c548851cc45b7b12d1707865`.
The workflow used MATLAB R2025b, pinned WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab`, and pinned Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Its
[runner](../../../matlab/ellipsoid_ode45_continuous_nonlinear_mass_control.m)
copied the published input and model, made the restoring, Froude–Krylov, and
wave-elevation blocks continuous, and used the [archived 0.05 s mass
control](../ellipsoid_ode45_continuous_nonlinear_mass_control/README.md).
Only the copied model's `ode45` MaxStep changed to 0.0125 or 0.00625 s
here. Output, nonlinear, and CIC clocks remained 0.05 s; duration remained
150 s. The 100 ns Transport Delay block remained present.

The mass control sets `adjMassFactor = A_zz / trace(A_xyz)` from the pinned
HDF5. It moves `A_zz` into the rigid inertial mass, leaving zero active
heave acceleration coefficient applied through the delay. Across all three
MaxSteps, the mass values were unchanged:

| Source mass identity | Regular | RegularCIC |
| --- | ---: | ---: |
| `adjMassFactor` | 0.655990137167053 | 0.783566018955714 |
| Original rigid mass (kg) | 133,376.067 | 133,376.067 |
| Heave added mass `A_zz` (kg) | 199,883.891 | 149,463.860 |
| Adjusted inertial mass (kg) | 333,259.957 | 282,839.927 |
| Residual active heave coefficient (kg) | 0 | 0 |
| Largest applied heave added-mass force (N) | 5.99e-22 | 3.76e-22 |

Each NPZ stores both refined steps on the same 3,001-point, 0.05 s output
grid: heave, speed, total and component source forces, acceleration, active
and inactive added-mass evidence, PTO stroke/speed/force/absorbed power,
incident wave, nine checked settings, and twelve mass values. The
[packer](../../../../tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control_step.py)
compared all 38 datasets of every exported HDF5 against the pinned source
HDF5 by value, including NaN locations. The saved wave agrees with both the
0.05 s control and paired delayed source within `1e-10` m. On every saved
state, the source mesh reconstructs restoring and viscous drag within
`7e-9` and `7e-11` N. PTO channels obey the 1.2 MN s/m damping and
absorbed-power laws. Inactive motion contributes less than `1e-6` N of
cross-axis added-mass force.

The following peak absolute gaps use **all 3,001 samples**, and work is
signed trapezoidal absorbed energy over 150 s. The independent Python run
uses the pinned HDF5/STL and published wave/PTO inputs, with no source
trajectory supplied to its dynamics. The paired delayed arms are the
[continuous-nonlinear 0.05 s base](../ellipsoid_ode45_restoring_timing/README.md)
and [refined source fixtures](../ellipsoid_ode45_restoring_step/README.md)
at the same MaxStep.

| Case | MaxStep (s) | No-delay versus Python: heave (mm) | Speed (mm/s) | PTO force (kN) | Power (kW) | Signed work |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Regular | 0.05 | 7.978 | 13.356 | 16.027 | 13.663 | +2.1080% |
| Regular | 0.0125 | 1.916 | 3.356 | 4.027 | 3.156 | +0.5352% |
| Regular | 0.00625 | 0.956 | 1.716 | 2.060 | 1.640 | +0.2688% |
| RegularCIC | 0.05 | 8.056 | 15.063 | 18.075 | 15.417 | +2.1690% |
| RegularCIC | 0.0125 | 2.073 | 3.699 | 4.439 | 3.410 | +0.5131% |
| RegularCIC | 0.00625 | 1.122 | 1.843 | 2.211 | 1.672 | +0.2618% |

| Case | MaxStep (s) | No-delay versus paired delayed source: heave (mm) | Speed (mm/s) | PTO force (kN) | Power (kW) | Signed work |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Regular | 0.05 | 4.661 | 7.538 | 9.046 | 7.762 | +1.2171% |
| Regular | 0.0125 | 1.100 | 1.853 | 2.223 | 1.764 | +0.2991% |
| Regular | 0.00625 | 0.537 | 0.951 | 1.141 | 0.916 | +0.1459% |
| RegularCIC | 0.05 | 3.880 | 7.074 | 8.489 | 7.217 | +1.0432% |
| RegularCIC | 0.0125 | 0.924 | 1.699 | 2.039 | 1.548 | +0.2300% |
| RegularCIC | 0.00625 | 0.447 | 0.850 | 1.020 | 0.778 | +0.1204% |

All motion, PTO, and work gaps contract over these three steps. The
0.0125-to-0.00625 no-delay change is still 0.977 mm heave, 1.784 mm/s
speed, 2.141 kN PTO force, 1.574 kW power, and 0.2650% signed work for
Regular; for RegularCIC it is 1.073 mm, 1.982 mm/s, 2.378 kN, 1.738 kW,
and 0.2500%. These finite-step observations do not establish an
asymptotic limit or isolate the 100 ns delay: `adjMassFactor` also shifts
inactive inertial masses and changes the solver's feedback representation.
The source restoring block still uses the original 133,376.067 kg
equilibrium mass.

At 0.00625 s the no-delay arms meet the inherited pointwise Python
ceilings (Regular: 3 mm, 4 mm/s, 5 kN, 4 kW; RegularCIC: 3.5 mm, 4.5 mm/s,
6 kN, 5 kW), but **both fail the 0.1% signed-work ceiling**. They are
post-screen diagnostics. They do not change the failed prospective source
step screen or the unchanged published-input parity gate.

The [portable PR
workflow](../../../../.github/workflows/ellipsoid-ode45-continuous-nonlinear-mass-control-step-portable.yml)
runs the [diagnostic test](../../../test_ellipsoid_ode45_continuous_nonlinear_mass_control_step.py)
without MATLAB. The two archives are 484,595 and 593,172 bytes, with
SHA-256 hashes `5d9eaefe9f305c0e5e7b41d811a232dd246ade78384e304601dd979ce079d4da`
and `692ffcb3ea207743ff7b77b7b0072a8ce3319bc4c382edf963c306c57cd19282`.
Repack and compare the downloaded artifacts:

```sh
gh run download 40131806849 -n ellipsoid-ode45-continuous-nonlinear-mass-control-Regular-maxstep-0.0125 -D .evidence/mass_control_step/Regular_0p0125
gh run download 40131806849 -n ellipsoid-ode45-continuous-nonlinear-mass-control-Regular-maxstep-0.00625 -D .evidence/mass_control_step/Regular_0p00625
gh run download 40131806849 -n ellipsoid-ode45-continuous-nonlinear-mass-control-RegularCIC-maxstep-0.0125 -D .evidence/mass_control_step/RegularCIC_0p0125
gh run download 40131806849 -n ellipsoid-ode45-continuous-nonlinear-mass-control-RegularCIC-maxstep-0.00625 -D .evidence/mass_control_step/RegularCIC_0p00625
.venv/bin/python tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control_step.py Regular .evidence/mass_control_step .evidence/regular-mass-step-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_continuous_nonlinear_mass_control_step/regular.npz
.venv/bin/python tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control_step.py RegularCIC .evidence/mass_control_step .evidence/regular-cic-mass-step-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_continuous_nonlinear_mass_control_step/regular_cic.npz
.venv/bin/python -m pytest -q -s tests/test_ellipsoid_ode45_continuous_nonlinear_mass_control_step.py
```
