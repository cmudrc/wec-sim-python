# Ellipsoid ode45 continuous nonlinear, no active heave delay control

This is a diagnostic for [issue #599](https://github.com/cmudrc/wec-sim-python/issues/599), separate from the unchanged published-input parity gate. The [R2025b source run 40130836521](https://github.com/cmudrc/wec-sim-python/actions/runs/40130836521) completed both 150 s wave cases at repository commit `87d27eb62b0d13d84247e113f22bae5d1fd60c03`. It pinned WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The [MATLAB script](../../../matlab/ellipsoid_ode45_continuous_nonlinear_mass_control.m) copies each application input and model into a disposable sibling directory. The original upstream inputs and models remain untouched.

The control retains `ode45`, 150 s duration, 50 s ramp, 0.05 s solver maximum/output/nonlinear/CIC clocks, the 4 m/6 s wave, mesh, hydro data, PTO, and the 100 ns Simulink Transport Delay. Like the [archived continuous-nonlinear source arm](../ellipsoid_ode45_restoring_timing/README.md), its three nonlinear blocks compile at a continuous sample period. The sole additional input change sets `adjMassFactor = A_zz / trace(A_xyz)` from the same pinned HDF5. In this heave-only model, WEC-Sim moves `A_zz` into the rigid inertial mass and leaves zero active heave acceleration coefficient entering the delay.

| Source mass identity | Regular | RegularCIC |
| --- | ---: | ---: |
| `adjMassFactor` | 0.655990137167053 | 0.783566018955714 |
| Original heave rigid mass | 133,376.067 kg | 133,376.067 kg |
| Heave added mass `A_zz` | 199,883.891 kg | 149,463.860 kg |
| Adjusted inertial mass `m + A_zz` | 333,259.957 kg | 282,839.927 kg |
| Residual active delayed coefficient | 0 kg | 0 kg |
| Maximum actual applied heave added-mass force | `5.99e-22` N | `3.69e-22` N |

All 3,001 samples share the same incident wave as the archived source arm and independently run Python within `1e-10` m. The active PTO speed equals heave speed and obeys the 1.2 MN s/m damping and absorbed-power laws. Inactive speeds and accelerations remain below `1e-9` in their units; the remaining off-diagonal heave mass couplings therefore apply less than `1e-6` N. The two regenerated HDF5 files agree in all 38 datasets with the [pinned source HDF5](../ellipsoid_ode45_convergence/hydroData/ellipsoid.h5).

These are peak differences over the entire 150 s between the matched mass control and the archived continuous-nonlinear source run. Signed PTO work rose by 1.217% and 1.043%, respectively, relative to that source run.

| Control minus archived source, peak magnitude | Regular | RegularCIC |
| --- | ---: | ---: |
| Heave | 4.661 mm | 3.880 mm |
| Speed | 7.538 mm/s | 7.074 mm/s |
| Logged total heave force | 302.729 kN | 173.441 kN |
| Restoring force | 9.776 kN | 9.215 kN |
| Excitation force | 8.303 kN | 8.465 kN |
| Radiation force | 0.884 kN | 0.608 kN |
| Heave acceleration | 0.03321 m/s² | 0.03476 m/s² |
| PTO force | 9.046 kN | 8.489 kN |
| Absorbed power | 7.762 kW | 7.217 kW |

The total-force change includes the shifted added-mass term; it is not a change in the pressure-force law at a fixed state. On every control state, the original-mass mesh reconstructs the source restoring and drag forces within `6e-9` and `7e-11` N, respectively. The source restoring block still uses the original `133,376.067` kg equilibrium mass. The change also shifts the inactive surge and sway inertial masses and changes the solver's feedback representation, so this control cannot allocate a numerical effect uniquely to the 100 ns delay.

The [portable Python test](../../../test_ellipsoid_ode45_continuous_nonlinear_mass_control.py) advances independent motion and PTO from the pinned HDF5/STL and published input. It compares all 3,001 source samples with that trajectory and with the archived continuous-nonlinear run:

| Case and source arm | Heave | Speed | PTO force | Absorbed power | Signed work versus Python |
| --- | ---: | ---: | ---: | ---: | ---: |
| Regular continuous nonlinear | 3.421 mm | 6.021 mm/s | 7.225 kN | 5.901 kW | +0.880% |
| Regular no active delay | 7.978 mm | 13.356 mm/s | 16.027 kN | 13.663 kW | +2.108% |
| RegularCIC continuous nonlinear | 4.307 mm | 8.196 mm/s | 9.835 kN | 8.199 kW | +1.114% |
| RegularCIC no active delay | 8.056 mm | 15.063 mm/s | 18.075 kN | 15.417 kW | +2.169% |

Removing the active delayed heave feedback **worsens every inherited motion, PTO, and work gap**. Both controls still fail the established 3/3.5 mm heave, 4/4.5 mm/s speed, 5/6 kN PTO-force, 4/5 kW power, and 0.1% work limits. This rules out the active delayed feedback as a sufficient explanation for the remaining continuous-nonlinear residual. It does not establish the underlying cause or promote either published input to physical parity.

`regular.npz` and `regular_cic.npz` retain time, active and inactive motion, source force terms, actually applied heave added force, PTO channels, wave, nine checked settings, and the twelve-value mass identity. Their SHA-256 hashes are `69018348739432556f08ad9e91bb32d855175b7b89e304174b8d9082309a98c7` and `364ef737b81a9ba3d2fd3b72ca82b0f3a95f68f42aa4aa16629ec6d6e0cfb435`. Regenerate and verify them with:

```sh
gh run download 40130836521 -n ellipsoid-ode45-continuous-nonlinear-mass-control-Regular -D .evidence/mass_control/Regular
gh run download 40130836521 -n ellipsoid-ode45-continuous-nonlinear-mass-control-RegularCIC -D .evidence/mass_control/RegularCIC
.venv/bin/python tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control.py Regular .evidence/mass_control/Regular/matlab-ellipsoid-ode45-continuous-nonlinear-mass-control/Regular .evidence/mass_control/regular-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_continuous_nonlinear_mass_control/regular.npz --source-hydro .evidence/mass_control/Regular/applications/Nonlinear_Hydro/hydroData/ellipsoid.h5
.venv/bin/python tools/pack_ellipsoid_ode45_continuous_nonlinear_mass_control.py RegularCIC .evidence/mass_control/RegularCIC/matlab-ellipsoid-ode45-continuous-nonlinear-mass-control/RegularCIC .evidence/mass_control/regular-cic-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_continuous_nonlinear_mass_control/regular_cic.npz --source-hydro .evidence/mass_control/RegularCIC/applications/Nonlinear_Hydro/hydroData/ellipsoid.h5
.venv/bin/python -m pytest -q tests/test_ellipsoid_ode45_continuous_nonlinear_mass_control.py
```
