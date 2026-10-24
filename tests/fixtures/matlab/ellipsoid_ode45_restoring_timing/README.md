# Ellipsoid ode45 nonlinear-force timing, full 150 s

This fixture is a diagnostic of [issue #599](https://github.com/cmudrc/wec-sim-python/issues/599).
The [R2025b source run 40129984760](https://github.com/cmudrc/wec-sim-python/actions/runs/40129984760)
completed both Regular and RegularCIC jobs at Python-repository commit
`9ed1e2ef2775f5bba02ace7433c3308772f74ddd` on
`codex/ellipsoid-force-timing-599`. Its pinned source checkouts were WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The source runner is
[`ellipsoid_ode45_restoring_timing_smoke.m`](../../../matlab/ellipsoid_ode45_restoring_timing_smoke.m),
called with `duration='full'`; it copies the application input and model for
each arm. Nothing is written to either upstream repository.

Both cases use the original 150 s, 4 m/6 s regular sea, 50 s ramp, ellipsoid
mesh and hydro data, 1.2 MN s/m PTO damper, `ode45`, 0.05 s solver maximum,
and 0.05 s physical/output/nonlinear/CIC clocks. The three arms differ only
in compiled sample time within the copied model:

| Arm | Nonlinear restoring | Nonlinear Froude–Krylov | Nonlinear wave elevation |
| --- | ---: | ---: | ---: |
| `baseline` | 0.05 s | 0.05 s | 0.05 s |
| `continuous_restoring` | continuous | 0.05 s | 0.05 s |
| `continuous_nonlinear` | continuous | continuous | continuous |

Each archive retains all 3,001 shared output samples for time, heave position
and speed, restoring and excitation force, acceleration, PTO stroke/speed/force
and absorbed power, wave elevation, plus all nine checked source settings for
each arm. The copied baseline matches the previously archived unchanged
published trace. The shared incident waves match across all arms within
`1e-10` m. The active PTO vertical speed matches body heave speed within
`1e-9` m/s; its force and power satisfy the source damping law. The 38 HDF5
datasets freshly generated in each source job agree in value with the
[archived source HDF5](../ellipsoid_ode45_convergence/hydroData/ellipsoid.h5),
although binary file hashes differ. The same archived source
[`elipsoid.stl`](../ellipsoid_ode45_convergence/geometry/elipsoid.stl) is used
by the portable Python test.

The source restoring log has the opposite sign to the mesh buoyancy-minus-weight
force reconstructed in Python. Over all 150 s, the baseline log follows the
**preceding** output state within `6e-9` N; each edited arm follows the
**current** state within `6e-9` N. Applying the opposite timing law leaves
over 51 kN of mismatch. This confirms that the copied-model intervention
changed the intended active force path. It does not change the original
published source input.

The [portable gate](../../../test_ellipsoid_ode45_restoring_timing_portable.py)
independently advances Python at the published 0.05 s step using only the
archived source-generated HDF5/STL and the configured wave/PTO inputs. It
checks all 3,001 samples and integrated signed absorbed work. These are the
measured peak gaps to that independent physical Python run:

| Case and source arm | Heave | Speed | PTO force | Absorbed power | Work |
| --- | ---: | ---: | ---: | ---: | ---: |
| Regular baseline | 18.940 mm | 28.341 mm/s | 34.009 kN | 26.270 kW | 5.819% |
| Regular continuous restoring | 5.156 mm | 5.860 mm/s | 7.032 kN | 7.216 kW | 1.433% |
| Regular continuous nonlinear | 3.421 mm | 6.021 mm/s | 7.225 kN | 5.901 kW | 0.880% |
| RegularCIC baseline | 19.930 mm | 30.042 mm/s | 36.050 kN | 28.771 kW | 6.194% |
| RegularCIC continuous restoring | 6.390 mm | 7.516 mm/s | 9.019 kN | 9.098 kW | 1.772% |
| RegularCIC continuous nonlinear | 4.307 mm | 8.196 mm/s | 9.835 kN | 8.199 kW | 1.114% |

Both edits improve every listed gap by more than half against the copied
baseline, but **both fail every inherited pointwise limit**. Regular's limits
are 3 mm heave/PTO stroke, 4 mm/s speed, 5 kN force, and 4 kW power;
RegularCIC's are 3.5 mm, 4.5 mm/s, 6 kN, and 5 kW. Each also misses the
unchanged 0.1% signed-work limit. The portable test treats the measured
improvement as a regression screen, not as a prospectively declared passing
threshold. The edited arms have no established source-step stability screen;
they cannot replace the separate derived fine-step reference. The unchanged
published-input gate stays red and both case rows remain `qualified_motion`.

`regular.npz` is 504,271 bytes (SHA-256
`0477f8433ff9dc4cb75913b36a726ab2a04478aa01137004f8fc3fbd85995a33`);
`regular_cic.npz` is 612,579 bytes (SHA-256
`f4d76114d8ba54c37a3256a19b2bec13a8a6ed754a37824c98cbb4626a96ccdd`).
Regenerate or compare each fixture from its downloaded workflow artifact:

```sh
gh run download 40129984760 -n ellipsoid-ode45-restoring-timing-Regular-full -D .evidence/force_timing_full/ellipsoid-ode45-restoring-timing-Regular-full
.venv/bin/python tools/pack_ellipsoid_ode45_restoring_timing.py Regular .evidence/force_timing_full/ellipsoid-ode45-restoring-timing-Regular-full/matlab-ellipsoid-ode45-restoring-timing/Regular .evidence/regular-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_restoring_timing/regular.npz --source-hydro .evidence/force_timing_full/ellipsoid-ode45-restoring-timing-Regular-full/applications/Nonlinear_Hydro/hydroData/ellipsoid.h5
gh run download 40129984760 -n ellipsoid-ode45-restoring-timing-RegularCIC-full -D .evidence/force_timing_full/ellipsoid-ode45-restoring-timing-RegularCIC-full
.venv/bin/python tools/pack_ellipsoid_ode45_restoring_timing.py RegularCIC .evidence/force_timing_full/ellipsoid-ode45-restoring-timing-RegularCIC-full/matlab-ellipsoid-ode45-restoring-timing/RegularCIC .evidence/regular-cic-repacked.npz --expect tests/fixtures/matlab/ellipsoid_ode45_restoring_timing/regular_cic.npz --source-hydro .evidence/force_timing_full/ellipsoid-ode45-restoring-timing-RegularCIC-full/applications/Nonlinear_Hydro/hydroData/ellipsoid.h5
.venv/bin/python -m pytest -q tests/test_ellipsoid_ode45_restoring_timing_portable.py
```
