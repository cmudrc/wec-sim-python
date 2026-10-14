# Portable MATLAB reference data

The published 100 s WaveBot adaptive controller input and trajectory are
documented in [wavebot_control/README.md](wavebot_control/README.md).
The published 600 s WaveBot impedance experiment is documented in
[wavebot_impedance/README.md](wavebot_impedance/README.md).
The two configured 400 s Sphere PTO comparisons are documented in
[sphere_pto_config/README.md](sphere_pto_config/README.md).
The public RM3 direct linear generator's published 20 s source window is
documented in [rm3_direct_generator/README.md](rm3_direct_generator/README.md).
The full 900 s published Sphere variable-mass comparison is documented in
[sphere_variable_mass/README.md](sphere_variable_mass/README.md).

## Wave-class outputs

`wave_reference.zip` contains 36 numerical CSVs and `matlab_release.txt` from
the successful [R2025b wave-source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38324036361)
(run commit `fcec46b15eab1bc3f8e461aedde6ac7ebf0c7212`). The workflow pins
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`; regenerate the
files with `tests/matlab/regular_wave_reference.m`. A separate local R2025b
Update 7 run produced numerical arrays exactly equal to all 36 archived CSVs.
The ZIP SHA-256 is
`d9ceb5e28daa572fd2bd05341c9782728df53e408089e457dc3bc7900fdc0668`.
Ordinary Python CI checks regular and finite-depth waves and six PM/JONSWAP
spectral, phase, power, and marker-output configurations against the archive.
Set `WEC_SIM_MATLAB_REFERENCE_DIR` to use fresh MATLAB outputs for the same
tests. These are wave-preprocessing checks, not published coupled-motion gates.

## Sphere latching

`sphere_latching.npz` retains all 40,001 samples of the published 400 s
`Controls/Latching` source time, body heave and speed, controller force, and
controller power. It comes from the successful [R2025b source sweep](https://github.com/cmudrc/wec-sim-python/actions/runs/38532934332)
(`matlab-Sphere_Latching-baseline`, run commit
`f8f8f14b1e8c70260f783df85f94b3176185d3b8`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
The body CSV SHA-256 is `7f5d5f4ae49b9df3346800dff1123b3467d06d3f5e2ed8c86f4377da512a8002`;
the controller CSV is `fb7407195929649aa41fad17940486cd2fcd9b0f6ab1dacacf7725b6850b4e48`.
Load both CSVs as float64, verify that their time columns match, and save body
columns 0, 3, and 9 as `time`, `position`, and `velocity`, and controller
columns 3 and 9 as `force` and `power` with `numpy.savez_compressed`.
The archive SHA-256 is `2fd3aec85a62068cd4c979528471c3597aa184ec3395165b1aa23820b4130f75`.
The source HDF5 SHA-256 is `98c847c4d5f106455f542f469c519f5d71fbb700b2f35b11c7f95dab21300116`;
all 46 datasets match bundled `examples/data/sphere.h5`, including
NaN positions. Ordinary CI uses the bundled HDF5; the live source workflow
continues to use its generated HDF5 and fresh CSVs. Both paths run the same
full-duration motion, latch-event, force-law, and absorbed-energy gates.

`oswec_spline_interpolation_probe.csv` is a six-sample numerical probe made
with MATLAB R2025b Update 7 from bundled `examples/data/oswec.h5` (SHA-256
`28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`).
Regenerate it with `tests/matlab/oswec_spline_interpolation_probe.m`. It uses
three between-knot frequencies, two BEM headings, fixed phases, MATLAB
`interp1(...,'spline')`, and a 0.2 s ramp. The portable test compares all six
excitation channels within `1e-7` N; the measured maximum error was
`3.8e-10` N. This is an interpolation check, not a published-case motion gate.

## Fixed hydrodynamic monopile

`monopile_hydro.h5` is the exact generated six-DOF HDF5 from the successful
[400 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37858193863)
(`matlab-MONOPILE_HYDRO-baseline`, commit
`94892a712e28c9243c6e7d60545800eb6136869d`). Its SHA-256 is
`927748f615f5e84e325e4d4ed96a01297dc5af40066f8c2c131e7e8d94930471`.
The workflow pins WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`monopile_hydro.npz` retains all 40,001 published samples of wave elevation,
both bodies' positions, velocities, and total forces, the hydro body's
excitation, and both fixed-joint reactions. It also keeps all 500 equal-energy
wave components. Arrays were loaded as float64 from the source CSVs without
rounding or interpolation; the archive SHA-256 is
`de4293b4e00a54a9bf72eba495b3a9af2273a9058810356dea5250fdd99a538d`.
Its source CSV SHA-256 hashes are:

| Source CSV suffix | SHA-256 |
| --- | --- |
| `components` | `df04a7be1b4516cdf7e686a4d528212e1cf3de525e592dd4ef69898d4882bc14` |
| `wave` | `d3ef951bde17a24fa47c22b7ce9bb34f90c1e34628520130c9a4e10c2c4381ab` |
| `monopile_body1` | `137d917ebd8b513d8d1b5471a5040df479d5d6e22c7b1a2e730c779dd8d7227f` |
| `monopile_body2` | `25aa640398864c7800b2a05e612aee9289161a1568e52bb4429e7987a5a152ad` |
| `constraint1` | `43520af1187b99bcce8c6af6ba7483f660ef826e002d1c3b496c0b1bb5aabce4` |
| `constraint2` | `07f081b62f722c82b3dd8046089b26d5473c9905c63bbf755e295ac48481fc06` |
| `body1_forces` | `f32e9da4f3fa2eb670dad15048df1f4679813cf5ea7a561fe2e97a7421afc98a` |
| `body2_forces` | `1c6675983839239d85436b6579d7c6356684d2a15faf7e848287b8ab56868ee4` |

To regenerate, load those CSVs with `numpy.loadtxt(..., delimiter=",")`.
Use wave columns 0–1, body columns 1–6, 7–12, and 13–18, hydro-body
excitation columns 19–24, and fixed-joint columns 19–24. Stack body and joint
arrays with body 1 before body 2; save them with `numpy.savez_compressed` under
the keys used by `test_monopile_hydro_portable.py`. The two `body*_forces`
files are not stored: the source test checks their stationary radiation,
added-mass, drag, damping, and acceleration channels are zero before packing.
The portable gate regenerates the sea from the pinned HDF5 and source phases,
then independently computes all 400 s of excitation and fixed-body loads.

## OSWEC nonlinear visualization motion

`oswec_nonlinear_viz_published_motion.npz` keeps all 1,201 samples of the
published 120 s, 0.1 s `ode4` flap motion and its principal force channels.
It was copied without rounding from the untouched R2025b
[published-case run](https://github.com/cmudrc/wec-sim-python/actions/runs/38259283403)
(`matlab-oswec-nonlinear-viz-pressure` artifact, run commit
`1083a1f29725785cbba0cb67f8e389a9cade16fc`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The source `source.mat`
hash is `f5c2944777f35741a48034e64fb4d3d7ec165466f6f2fc76f7be1f38383ad466`.
The archive contains `time`, `pose`, `velocity`, `acceleration`,
`forceExcitation`, `forceRadiationDamping`, `forceAddedMass`, and
`forceRestoring` for the six flap degrees of freedom. Its SHA-256 is
`63f3669b192ab271a92f5ea9c4650d1d183ef68f0abc4a6b0faedcd2b4948e4f`.
The refined-step source job gates its 0.1 s run against these published
channels before evaluating step sensitivity or independent Python motion.
The motion archive does not contain hydrodynamic coefficients. The pinned
[refinement artifact](https://github.com/cmudrc/wec-sim-python/actions/runs/38093454463)
supplies its matching `oswec.h5` (SHA-256
`72664d60037698350513d6481d4114db79fe2cf8984247a13ff16e980c52d31f`).
That HDF5 has 36 excitation headings; `examples/data/oswec.h5` has 10 and
also differs in added mass. Use the source-matched HDF5 for this motion
comparison rather than substituting the example file.

## RM3 regular wave

`rm3_regular_wave.npz` contains selected motion and PTO channels at every
sample of the 400 s published RM3 regular-wave baseline. It was extracted
without rounding from the successful
[pinned R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/38212367776)
(`matlab-RM3-baseline` artifact, run commit `dbcd62bfeba7bcc857f699522d17c06add834c39`).
The workflow checks out WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab`. Inputs are a 2.5 m,
8 s regular wave, 0.1 s output step, 100 s ramp, and a 1.2 MN s/m linear PTO.

The archive contains `time`, `body_position`, `body_velocity`,
`body_heave_excitation_force`, `pto_stroke`, `pto_speed`, `pto_force`, and
`pto_power`. Body arrays are ordered float/spar, then time; position and
velocity add surge/sway/heave/roll/pitch/yaw as the last axis. These are
source CSV columns 0, 1–6, 7–12, and 21 for bodies, and columns 3, 9, 15,
and 21 for the PTO.
The source files in the artifact have these SHA-256 hashes:

| Source file | SHA-256 |
| --- | --- |
| `RM3_RM3_body1.csv` | `f3f31922d9fac06b13b1ad6393443256084d30c86ffe1788eca4a97e3273ff50` |
| `RM3_RM3_body2.csv` | `60594f5c7cebf17e3dd8347210abab4b1324b3e3ba50b9cc8631c7b66d88ecaf` |
| `RM3_RM3_pto1.csv` | `d67f8833edb5cd67bdba5a48b76a6a440cd7eed68dce4800a108d478a014f600` |

To regenerate, load each CSV with `numpy.loadtxt(..., delimiter=",")`, check
that all three time columns match, and save those columns using
`numpy.savez_compressed` under the names and ordering above. The archive SHA-256
is `bc8db627eb710b002e6ed91b0a352999ea62a56c09baf5e859677fa1f1442ba5`.

The always-on test uses the bundled `examples/data/rm3.h5` (SHA-256
`3bdfdcfa3dcd2db1851ab2026eb6b4fdc565ad86be8503e09b42231feee429fb`).
That HDF5 is not byte-identical to the current WEC-Sim example; its regular-wave
hydrodynamic data still reproduce this source trajectory within the stated
per-channel gates. The live paired tests continue to use pinned MATLAB inputs.

## OSWEC Wave Markers

`oswec_wave_markers.npz` retains every 0.1 s sample of the published 400 s
OSWEC regular-wave marker, body-motion, and PTO histories. It was copied
without rounding from the successful
[R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068)
(`matlab-WAVE_MARKERS_OSWEC-baseline` artifact, run commit
`9872f98ccd65aa57e49f0d7dbd5db8919d87773c`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The source marker CSV hash
is `90aab5317195b8fcb656089fdd42e2de0ccae2d96b1fb948b5b2a9916acef8c3`;
the archive SHA-256 is
`16709146ff26f0d5b971f16a48713af8dd58ea325870802732cd1d4db148c2b1`.

The archive stores time, 25 world-XY marker locations and elevations, origin
elevation, the seven source wave parameters, both bodies' six position and
velocity channels, and hinge PTO stroke, speed, force, and source-signed power.
For each body, position and velocity come from source CSV columns 1–6 and
7–12. PTO stroke, speed, force, and power come from columns 5, 11, 17, and
23. The source body CSV SHA-256 hashes are
`3e6363059307858a8795e2e9390430e8fa753f0dc15ba6eb362c909296cc8b2f`
(flap) and
`7c4dec9675b5fd85007259728902aeaf37c463c266ccfdbccd043e77bb38e37d`
(base); the PTO CSV hash is
`b297702932cf6147efccb64b462036064d1b26a8e39e9795fd593b3ad62dd4c7`.
The always-on test uses `examples/data/oswec.h5` (SHA-256
`28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`).
This bundled HDF5 is not byte-identical to the source job's generated HDF5
(SHA-256 `79d9ab232958f3b3f9ec1bcdcd9a56c9fe8eb54f22e9de2714403860b8a51b32`),
but its full-duration motion passes the same per-channel source limits. The
external paired test still checks the exact generated HDF5; this portable
gate keeps the measured published trajectory in ordinary Python CI.

## GBM Barge

`gbm_barge.h5` is the exact generated ten-DOF HDF5 from the successful
[400 s R2025b baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/38282150325)
(`matlab-GBM_BARGE-baseline` artifact, run commit
`7c403d3e2e07b9940e4a2a280d375000b8985fea`). Its SHA-256 is
`0d98003d54f007f1c50966bb10df723f5bbbf22aa70579360e4ef991a2ceff68`.
The workflow pins WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`gbm_barge.npz` keeps all 8,001 source samples of rigid body time, pose,
speed, and acceleration (CSV columns 0–12 and 43–48), the 40-channel
`Flex_out` signal, and the wave, mass, stiffness, damping, and hydrodynamic
matrices used by the existing full-duration gate. Arrays were
loaded from the source CSV and MAT files as float64 without interpolation or
rounding. The `Flex_out` MAT SHA-256 is
`2487a21c2e8760bea072f39ec9205c7992927593c7c1704714f907c9fbd59709`,
the rigid body CSV SHA-256 is
`666e42a9c6ce0e01e12dc6a35b20d542f7531682a2882b41b29b209f98721ecf`,
and the compact archive SHA-256 is
`e9c8cddf5f10c9d60ddf26348082a36cd436fbdaf59decac4e4a2870797df1c9`.
Ordinary Python CI uses these files when live MATLAB output is unavailable;
the same test uses fresh generated files when both source directories are set.

## OWC OrificeModel physical prefix

`owc_orifice.h5` is the exact seven-DOF BEMIO output from the successful
[R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37799665808)
(`matlab-OWC_ORIFICE-baseline`, commit
`5dbf0d5fc228bcbd39275ec8c267ef0e43289bb2`). Its SHA-256 is
`fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`.
The workflow pins WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`owc_orifice_prefix.npz` keeps the first 1,201 source samples (0–6 s) of
rigid-body output, flexible-mode signals, and incident sea, plus all 500 PM
component rows, the five published orifice parameters, and five unrounded
internal orifice evaluations. Those evaluations cover both flow signs and
subthreshold/superthreshold Mach flags; their original row indices are stored.
The rigid-body CSV SHA-256 is
`8153820cc0de55b347ca0f7330eea49d1e5c9e71176c8145898c5687b97219b2`,
the flexible-mode MAT SHA-256 is
`825ff26839df1fccb5c48533ec71511671240d5693450cf0cb99df1518c6e20c`,
and the full internal orifice CSV SHA-256 is
`22314e78887134feebc28f2a54e1e5157d39769ce534229367f44dc22cd43096`.
The compact archive SHA-256 is
`d6106f4b14abfed868665f478f2b658346f535f501c5ec533e356620bc79dca7`.
The ordinary-CI Python run receives only the pinned HDF5 and saved wave
phases. Its motion is compared after the solve against the same predeclared
6 s limits as the live source test. This pre-threshold prefix does not gate
the published 130 s physical coupled trajectory.

## MBARI Cable

`mbari_cable.npz` keeps all 12,001 samples of the published 120 s three-body
MBARI Cable motion, incident wave, and cable displacement, speed, and axial
force. It comes from the successful
[pinned R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37859525482)
(`matlab-CABLE_SOURCE-baseline` artifact, run commit
`fc753cf65892553cf3363b496ae7111444b096e5`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The matching generated
`mbari.h5` was copied unchanged to `examples/data/mbari.h5` (SHA-256
`29da78db89775f682ba1b50b613f594e4a961195529157967493b3ffb8794422`).

The archive retains body CSV columns 1–6 and 7–12 as position and velocity,
wave columns 0–1 as time and elevation, and cable columns 3, 9, and 21 as
displacement, speed, and force. All columns were loaded as float64 without
rounding. Source SHA-256 hashes are:

| Source file | SHA-256 |
| --- | --- |
| `CABLE_SOURCE_Cable_body1.csv` | `ec8bd074a299a0a3bc0e792432495a5b0c26abbc92aa974198853ab21a86f3b0` |
| `CABLE_SOURCE_Cable_body2.csv` | `42c2b464b93099b57b88a00d0a3f45d2b7f1b8015706db81483ebd37ed1c6e4a` |
| `CABLE_SOURCE_Cable_body3.csv` | `61f95b211e3d34bbbef0543c6b6f3ec12e0eb110e6ac151d6395a5b625d82561` |
| `CABLE_SOURCE_cable1.csv` | `aca36883789e40e53476123f10b45b3db91e768f24e657f66add7d5bd9e27859` |
| `CABLE_SOURCE_wave.csv` | `de635cbf9f60210cdaeab3ce89782e3b657df85a3a682dac809f0af2f17bc925` |

The compact archive SHA-256 is
`dfbe2cd5d689b32c436d243e0772d6298632bc6221bfe45c40a83a0d6d882e7e`.
The always-on test independently advances the complete published input and
applies the established body, wave, cable, force peak/RMS/impulse, and active
fraction gates. This guards the published planar case, including its brief
step-sensitive cable snap; it does not establish arbitrary six-DOF cables.

## WaveStar base case

`wavestar_base.npz` retains all 14,121 samples of the published 141.2 s
`WECCCOMP/WECCCOMP` float, arm, PTO, excitation, radiation, wave, and
500-component sea. It comes from the successful
[pinned R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37866626954)
(`matlab-WECCCOMP_SOURCE-baseline` artifact, run commit
`d4dfb38c230a31e07d4525671e758e49f4365bc9`). The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Its generated HDF5 is
copied unchanged to `examples/data/wavestar.h5` (SHA-256
`cf61dfada9cf4d177a02dbad8c4f4eaf7c4ac022fbf188faca9ccfb62a6a0faf`).

The archive retains the source component columns, body positions (columns
1–6) and velocities (7–12), float excitation (19–24), PTO stroke and speed
(3 and 9), radiation force (1–6), and wave time/elevation (0–1) as float64
without rounding. Source SHA-256 hashes are:

| Source file | SHA-256 |
| --- | --- |
| `WECCCOMP_SOURCE_components.csv` | `80bc9515aa5abcc8de88d01a74737ea9f059f4743cfe18fddb22d2cdd8630f9e` |
| `WECCCOMP_SOURCE_WECCCOMP_body1.csv` | `5067a2119405ce2cecc5a09ee6ca5ddd01867c4f03cec43740a45f6c4b408c82` |
| `WECCCOMP_SOURCE_WECCCOMP_body2.csv` | `4439a0a875fd68617efa9d5da41029fe1d15fa471a80a49d72258f5497217ed3` |
| `WECCCOMP_SOURCE_WECCCOMP_pto1.csv` | `0e0c678eb833b125d3f252a8a382691ac978f7c6559a06ce93690c7d28feb4f6` |
| `WECCCOMP_SOURCE_wave.csv` | `4c30e47d607abc1b680d26f1c61e8edc1b539b0efb52045eeecafd7583e96e38` |
| `WECCCOMP_SOURCE_body1_forces.csv` | `c23a085c7b444c90f6f8f5ab574de3bbe94807ad9fa2f040466280f51fd481d5` |

The compact archive SHA-256 is
`46bf4853f3bb8a5a5572adde2f2871259ee87fdc3f200938f3070e2b26c18f95`.
Ordinary Python CI independently integrates the full published base case from
the HDF5 and sea phases, with the established motion, PTO, and force gates.
It does not gate the separate fault or NMPC variants.

## WaveStar fault case

`wavestar_fault.h5` is the exact BEMIO-generated HDF5 from the successful
[141.2 s R2025b fault run](https://github.com/cmudrc/wec-sim-python/actions/runs/37871185462)
(`matlab-WECCCOMP_FAULT_SOURCE-baseline` artifact, run commit
`be6678cbd9f029c043dbfc61738534679b60dc21`). Its SHA-256 is
`a753474c44592684ba6eadac4d01a0e08ef0714f5e58a67383cb0dee9f70914e`.
The workflow pins WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.

`wavestar_fault.npz` retains all 14,121 saved output samples and all 141,201
sensor samples. It stores the 500 JONSWAP components, the measured-minus-true
sensor disturbance and dropout flags, float and arm surge/heave/pitch motion,
six float excitation and radiation channels, PTO stroke/speed/force, and wave
elevation. The source CSV columns were loaded as float64 without rounding;
the Boolean dropout is true exactly where the measured position is zero.
The original true and measured sensor CSV SHA-256 hashes are
`e72e091958d164a75dbf07c2d9c6833667700c1f3d2760293cc4fc683d7c2139`
and `18bd91ac4bb1e986c8c375940fc40cd2b9bb7809700d9c5910eb0df304ff6203`.
The float body CSV hash is
`672e9f7db671a3f737970eea4ff3ed9a6c4ef13185290f1f8e8e0766edce106b`;
the compact archive SHA-256 is
`9102b671cf31bfa41b6d91789d031beb1b34c95e942f2df62551a0f5ab1ff8d8`.
The always-on test independently advances the public fault solver from these
exogenous sensor and sea inputs. Saved MATLAB motion and force histories are
used only as expected outputs; the live-artifact mode remains available.

## RM3 MoorDyn visualization source-step audit

`rm3_moordyn_source_step.npz` retains the same-seed 80 s MATLAB source
trajectories from the [published 0.01 s run](https://github.com/cmudrc/wec-sim-python/actions/runs/38403345509),
the [MaxStep-only 0.0003125 s run](https://github.com/cmudrc/wec-sim-python/actions/runs/38117328859),
and the [fully refined 0.0003125 s run](https://github.com/cmudrc/wec-sim-python/actions/runs/38108596566),
in that order. Their repository commits are `f3f12b342d9f0005a1174013b4e34da34b16f8be`,
`865d6f8670e3670dded3e6cb38df7ad65216ebfc`, and
`05a0986d3597c9f67aac5e0733ca96e2fe7ee4ea`. All workflows pin WEC-Sim
`0753b2e47f2457c078751dcfe5d251d1767b80ab`, Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, MoorDyn
`34c0703c9a5b74fd2cc5bb52caaadaaa356d675e`, and MATLAB R2025b.

Each source `phase.csv` supplies all 1,000 phases. `wave.csv`, `body1.csv`,
`body2.csv`, and `mooring.csv` supply elevation, position columns 1–6, and
mooring-force columns 13–18. The two 0.01 s exports are sampled every tenth
row to align with the published 801-row, 0.1 s export. No source data are
interpolated or rounded. The archive SHA-256 is
`4f7f747465470caf1693b7a438e2bb40bd0835676c571a0230401bd80fba7db2`.
The portable test gates source input identity and the channel-specific effect
of ode45 MaxStep on the MATLAB trajectories. It is a source-numerics diagnostic,
not an independent Python published-step motion gate.
The same archived phases also drive an independent Python reconstruction of
the published 80 s incident sea in ordinary CI, without using saved motion.

## RM3 End Stops full-duration trajectory

`rm3_end_stops_full.npz` contains all 400 s of the pinned R2025b
`RM3_END_STOPS_FULL_FINE` and `RM3_END_STOPS_FULL_FINER` jobs in
[run 37698235851](https://github.com/cmudrc/wec-sim-python/actions/runs/37698235851)
(repository commit `ae9a0a8e0959b95a9489b61ea9c6f161ae184e7e`). The workflow pins
WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The finer job's generated
`rm3.h5` is copied unchanged to `examples/data/rm3_end_stops.h5` (SHA-256
`e76bff162db8ab1f93747abf35142e733d9e0c43cb9751aee5bb9e5636b87176`).
Its fitted radiation data differ from `examples/data/rm3.h5`, so that general
example file cannot stand in for this source run.

The archive retains the original run/stop settings, all 16,001/32,001 time
samples, each body's surge/heave/pitch position and speed (source columns
1/3/5 and 7/9/11), and PTO stroke/speed/force (columns 3/9/15) as float64.
It neither interpolates nor rounds source values. Archive SHA-256:
`4d21c43c70fc7c213136c73e1450c0480c0739c76d8c2c6910c3ed51f4c858e3`.
Ordinary Python CI now gates both the two-step MATLAB source stability and an
independent 400 s public `WEC.run` solve using only the generated HDF5, wave
settings, and configured stop law. This is the refined-step case; the
published 0.1 s source trajectory remains step-sensitive and qualified.

## RM3 end-stop force transition

`rm3_end_stops_force_samples.csv` retains six rows from the successful
[full-duration 0.0125 s R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37698235851)
(`matlab-RM3_END_STOPS_FULL_FINER-baseline` artifact, run commit
`ae9a0a8e0959b95a9489b61ea9c6f161ae184e7e`). The workflow pins core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Source rows
`0, 5166, 8627, 8948, 22725, 28451` retain CSV columns `0, 3, 9, 15`
as time, stroke, speed, and total PTO force. The source CSV SHA-256 is
`898147eeac1bc62746dfb84a71dec9d640bbcafb9da63129ef9040467398f7d0`;
the compact fixture SHA-256 is
`a409f621d21825366dc1f11b13c12fc9c88f07beb6f86a4670a709253c7f3c20`.
Two selected rows are inside the effective 0.1 mm spring transition and
distinguish the source's smooth force from a sharp stop.

## RM3 PTO extension source-step diagnostic

`rm3_pto_extension_step.npz` retains every 30 s float and spar free-decay
sample at 0.1, 0.05, 0.025, and 0.0125 s from the successful
[R2025b step audit](https://github.com/cmudrc/wec-sim-python/actions/runs/38324343546)
(commit `58ade4aa4ba8c4d94827696b6d722a42a5cf91e5`). Each nine-column CSV
is copied as float64 without interpolation or rounding. The archive SHA-256 is
`77a28920004ec030b2f34f5047d81848bc9aa8ba09f3d49fd7cefa8564a2060b`.

The portable diagnostic uses `examples/data/rm3_end_stops.h5`, also generated
from the pinned common RM3 BEM input. Its active infinite-frequency added
mass, radiation impulse response, restoring, displaced volume, and CG datasets
match the [archived PTO-extension baseline HDF5](https://github.com/cmudrc/wec-sim-python/actions/runs/38282150325)
exactly; their fitted state-space groups differ and are unused here. The test
checks source step sensitivity, public implicit-solver self-refinement, and
the opt-in delayed-mass numerical comparison. This remains a diagnostic; it
does not promote either published 0.1 s case.

`rm3_pto_extension_00625.npz` adds both complete 30 s trajectories at
0.00625 s from local MATLAB R2025b Update 7 with Simulink, Simscape, and
Simscape Multibody. The source checkouts remain at WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. For each case, a disposable
copy of the pinned `RM3_PTO_Extension` input changed only `simu.dt` from
0.1 to 0.00625 s. Its common-RM3 HDF5 was generated by that checkout's
`bemio.m` (SHA-256
`b8dd385f5db2a53fa3ffb971b105ea4fe8a2685ffd00bab0ae6050f366f69768`).
Before refinement, local 0.0125 s replays of both cases matched the saved
[R2025b audit](https://github.com/cmudrc/wec-sim-python/actions/runs/38381113784)
to at most `1.1e-13` m or m/s in active motion. The new 4,801-row CSVs use
the same nine columns as the earlier archive: time, float heave/speed,
spar heave/speed, and PTO stroke/speed/force/power. Their SHA-256 hashes are
`9483629afaee1e87c401b2dca76cfcff288ddd7a67401f2dc115d8736ac96b12`
(float) and
`7e5104cb01b277867f455c083384cde6fb6691e857f238a2e083394973637f4e`
(spar). Load them with `numpy.loadtxt(..., delimiter=",")` and save them
without rounding as `float_0p00625` and `spar_0p00625`; the compact archive
SHA-256 is
`c2428bd6488b5c90834d9a58dd068a7e7d7ae804688281c0fff1e5ef000fd4a3`.
The portable test compares the new source step change and independent public
Python motion. The published-step cases retain their qualified status.

`rm3_pto_extension_003125.npz` adds both 9,601-row, 30 s source trajectories
from a second local R2025b Update 7 refinement. Disposable copies of the
same pinned float and spar inputs changed only `simu.dt` from 0.1 to
0.003125 s, and reused the generated common-RM3 HDF5 with SHA-256
`b8dd385f5db2a53fa3ffb971b105ea4fe8a2685ffd00bab0ae6050f366f69768`.
The float/spar CSV SHA-256 hashes are
`2a3f9c92aec9bebccbeabc14cc7fba557cddefef93a1801276fb3fe95d5341f4`
and `7fb93a6366764409098e71049db797bcc4ee191d67ea8cf5add82d911dbb47e5`;
the compressed fixture SHA-256 is
`94ca53190aad68a0bd65bdefa9eca372e7873cdf5dfcc91a7339309fede16589`.
The always-on test uses the previously declared decreasing source-step and
extrapolation-inconsistency screens, then checks independent physical Python
motion against the latest formal extrapolation. This is a convergence
diagnostic, not a published-step parity gate.
The next two local source halvings and their independent physical comparison
are documented in
[rm3_pto_extension_00078125/README.md](rm3_pto_extension_00078125/README.md).

## Sphere reactive PI controller

`sphere_reactive_pi.npz` keeps every sample of the published 200 s Sphere
`Controls/Reactive (PI)` heave, speed, controller force, and controller power.
It comes from the successful [pinned R2025b baseline
run](https://github.com/cmudrc/wec-sim-python/actions/runs/38497841030)
(`matlab-Sphere_Reactive_PI-baseline` artifact, run commit
`397f48a71cf017b396a0e98e79ad63340749f778`). That workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The generated source
HDF5 is bundled at `sphere_pto_config/sphere.h5`. The archive stores CSV
columns 0, 3, and 9 from the body and columns 3 and 9 from the controller
as float64 without rounding. The source body and controller CSV SHA-256 hashes
are `34a9d0e44006f8c42ed8a866b20bf19f707f017ccc36865d0d0a7484072fcf15`
and `545c9cefdc26ca25e00e353ef42d1d349d9ab9ea210aa15d11fb1cb0bc9d195a`;
the compact archive SHA-256 is
`73bf7245e0b536afae5f64f515da0858f31520a220d32f3400ef23d0dffaa471`.
The test retains its 200 s motion, force, power, and MATLAB gain-law gates in
ordinary CI. This strengthens an already counted Sphere case.

## Sphere free decay

`sphere_free_decay.npz` retains every sample of the five published Sphere
`noWaveCIC` free-decay cases, ordered `0m`, `1m`, `1m-ME`, `3m`, `5m`.
The successful [pinned R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/38196174533)
(`matlab-Sphere-baseline` artifact, run commit
`33422a97e5e1273c2828b923d2f96d7c37f31caf`) used WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Its BEMIO-generated
`sphere.h5` is copied to `examples/data/sphere.h5` without modification
(SHA-256 `7287c5b02c74616e385f0d582df6a4679d7a4ceb93d8c663833a944838598e06`).

The archive stores `case_names`, a common `time`, and arrays ordered by case,
time, then WEC-Sim degree of freedom for `body_position` and `body_velocity`.
`body_total_heave_force` has case and time axes. Body CSV columns 0, 1–6,
7–12, and 15 were copied as float64 without rounding. Source hashes are:

| Source file | SHA-256 |
| --- | --- |
| `Sphere_0m_body1.csv` | `dd2f2f5d793c769df4ffe9610d4ad5773bf7d4780dab1b224babe88e10ba24ed` |
| `Sphere_1m_body1.csv` | `8eeacfa184077705ec88b8a87f35c32eff097c2a0ee344bca1216a0406ebf427` |
| `Sphere_1m-ME_body1.csv` | `8eeacfa184077705ec88b8a87f35c32eff097c2a0ee344bca1216a0406ebf427` |
| `Sphere_3m_body1.csv` | `4205f502b2c4c0fd14d14c1ba2c01a35bab5331ede5f1d36cb8cad9eb968e44b` |
| `Sphere_5m_body1.csv` | `6d66112742795e9d57009e08b1e88c5f68040669d5fdd589f7d962840d438876` |

The archive SHA-256 is
`95057bc8d58b1e36243eefb688a184774a557dc63542128efff3eac732095154`.
The two 1 m source trajectories are identical: the published Morison element
has surge-only coefficients and does not change heave free decay.

## Sphere moving-Morison heave

`sphere_moving_morison_heave.npz` keeps all 4,001 samples of the derived 40 s
`Free_Decay/1m-ME` case with its body-local element changed to axial heave
coefficients (`Cd=Ca=1`, area 100 m², volume 20 m³). This is the successful
[pinned R2025b MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37770321922)
(`matlab-SPHERE_MOVING_MORISON-baseline` artifact, run commit
`2c707cf2cd45e144db803a8862485d4845c6d992`). Its workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The bundled
`examples/data/sphere.h5` is the matching BEMIO output (SHA-256
`7287c5b02c74616e385f0d582df6a4679d7a4ceb93d8c663833a944838598e06`).

The archive stores source CSV columns 0, 3, and 9 as `time`, `position_heave`,
and `velocity_heave`; `physical_force_heave` is the negative of column 27,
because the MATLAB log uses the opposite Morison force sign. Values were
loaded as float64 without rounding. The source CSV SHA-256 is
`16788d9ae5d8d3413affadc722f2dc405f6f9cb76cecefa14ab1849a544ae0d8`;
the compact archive SHA-256 is
`3a35601b7b34103880082962ec828d64c82f26b710e32dda78e185f524141b7f`.
The always-on test independently integrates the complete motion and gates
heave, speed, force after startup, and the full force-impulse difference.
This derived heave-only experiment supplements the published Sphere cases;
it does not make a new published-case badge claim.

## OSWEC variable-yaw 160 s source step diagnostic

`oswec_yaw_160_source.npz` preserves the full 250 s `dt=0.005` MATLAB
trajectory and wave, plus every sample of the local 160 s `dt=0.0025`
trajectory and wave. The source is WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e` under local MATLAB
R2025b Update 7. The existing coarse/fine inputs came from the pinned
[250 s source audit](https://github.com/cmudrc/wec-sim-python/actions/runs/38127365821).
The committed [MATLAB runner](../../matlab/oswec_variable_yaw_160s_step_audit.m)
sets only `simu.endTime=160` and `simu.dt=0.0025` on the published irregular
input and uses the selected 0.25° heading bank. Its generated OSWEC HDF5
SHA-256 is `6b6a2837d16808382c4eeeaef5abe5cef5aab5a60cb8dd5e0124b9bac51b17f9`.

The archive retains unrounded float64 CSV values. SHA-256 hashes: source
`fine_trajectory.csv` `2276cace8ac5254002ff9e939ff3db51eeba417b1b6aa5f64003f9219802039d`,
`quarter_trajectory.csv` `375f9fb517f386257c5d538960630497111f4251dcca03c918c1e19f0dd05d2c`,
`quarter_wave.csv` `77188908365176bc5992a79df125291cc32d23ef8682e268ebcc6`,
and compact archive `dbf479b062288c19eafcd0d2f598e0838580f89774d52fce174131ca94ef3bcc`.
`oswec_yaw_160_metrics.json` contains the predeclared 120/140/160 s source
self-change and independent Python/source errors; its SHA-256 is
`aa44a6f24ddef3a8d8db30ebe517d7eedfb4b4505e04c43966efa025b31702a1`.
The unchanged limits and 160 s duration were recorded before the local run in
`oswec_yaw_160_gate_spec.json` (SHA-256
`a50186a2ce9b7ce9a611da5c95e9a2eb312991368b7255c3ef1b486311dc94bb`).
The Python comparator used the physical-default `WEC.run` path at repo
commit `547dd5ed` and the same generated HDF5; its 160 s comparison passed
the unchanged Python/source limits, while MATLAB's own step stability failed.
The always-on diagnostic recomputes the source side of that boundary from
the saved trajectories. This is not a published-case parity gate.

## OSWEC directional irregular wave

`oswec_irregular_wave.npz` keeps all 4,001 samples of the published 400 s
OSWEC PM sea, both body motions, the flap's six excitation forces, and the
hinge PTO. The source is the successful [pinned R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068)
(`matlab-WAVE_MARKERS_OSWEC-baseline` artifact, run commit
`9872f98ccd65aa57e49f0d7dbd5db8919d87773c`). It checked out WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab`. The published input
has three incident directions and no phase seed, so the archive retains its
500-by-3 phase realization. `examples/data/oswec.h5` is an unchanged copy of
that pinned core revision's `examples/OSWEC/hydroData/oswec.h5` (SHA-256
`28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`).

The archive stores `time`, `wave_elevation`, `wave_components` (frequency,
amplitude, width, then three phases), `wave_directions` (direction and
spreading), `body_position` and `body_velocity` (flap/base, time, six degrees
of freedom), `flap_excitation_force`, and PTO stroke, speed, force, and power.
The body CSV uses columns 0, 1–6, 7–12, and 19–24; the PTO uses columns 0,
5, 11, 17, and 23. All source CSV values were loaded as float64 without
rounding. Their SHA-256 hashes are:

| Source file | SHA-256 |
| --- | --- |
| `OSWEC_OSWEC_body1.csv` | `9b005121c637987aff76915f1b0295200970c7e3bec3711956e9cc25b8570edf` |
| `OSWEC_OSWEC_body2.csv` | `cc8d525d09cc876ae50766f6d1067face24f9a10a823afa19a06d702e99eb33b` |
| `OSWEC_OSWEC_pto1.csv` | `3a28e1d7dbcff757fb980f35ce828f895fb0b9f1481ca0f22ceb90e149775766` |
| `OSWEC_wave_components.csv` | `29a5ce24a901c819953da8ed365bc4c22e8b2beb66f2a8529124db6c19191b00` |
| `OSWEC_wave_directions.csv` | `fb581ce2fe25dcc31394e43c243ff2281627613fe0df5a543f712d4e9343cebf` |
| `OSWEC_wave_elevation.csv` | `35aae2f59b1b66aa2ab6107bca0eca712f2ffbb8c58bce208d2f9d9fcbfdc62b` |

The archive SHA-256 is
`070365a8fd15dddab8b4fe89a6000eb8572f291ee6f2068d3520e13612424b01`.
