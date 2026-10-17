# OSWEC hinge PTO pretension

`source.npz` retains the full 4,001-sample, 400 s float, fixed base, and PTO
histories from the [pinned R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/38201051995), plus its no-PTO flap control and four PTO parameters. The source used WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Each array was read from the run's `matlab-OSWEC_Nonhydro-baseline` artifact as float64 and compressed without rounding. SHA-256: `1d8bfaa3dbe807aa53e6684c4cab5f28ed09ca35bcb4c94b007c6a6223a3b897`.

The ordinary Python gate uses the already bundled
`../oswec_multi_wave/oswec.h5`. Its 68 datasets match the archived run's HDF5
except four impulse-response arrays, whose maximum element difference is
`4.01e-11`. Both public Python hinge paths pass the existing full-duration
position, speed, PTO torque, power, and spring-effect limits with either HDF5.
When fresh MATLAB CSVs are supplied, the same tests compare directly to them.

This is a derived PTO configuration. The source rotational mask interprets its
equilibrium offset in degrees; the Python test converts that value explicitly
and does not claim that the same literal pretension has the same units in both
APIs. This fixture adds a durable configuration gate, not a published-case
badge.
