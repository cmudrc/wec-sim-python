# OSWEC regular-wave solver generalization

This is a **derived** 60 s input, not an additional published case. It copies
the pinned Applications `Wave_Markers/OSWEC` model and changes only wave
height from 2.5 to 1.3 m, period from 8 to 6.5 s, ramp from 100 to 15 s,
duration from 400 to 60 s, and disables post-run explorer display. The
published `ode4` step remains 0.1 s. The HDF5 SHA-256 is
`28aa028d404a88f877924bc435e7ad3a4c50e26f9706c1f39ef7324e5ec60750`.

- WEC-Sim source: `0753b2e47f2457c078751dcfe5d251d1767b80ab`
- Applications source: `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`
- MATLAB: local R2025b on Mali, WEC-Sim 7.1.0
- `source.npz` SHA-256: `e1350e37d974f9d1f6ffac4c0e82f78104fda960a871367e38d16b7434cde552`

Run `oswec_regular_derived_60s.m` with the two pinned repository paths,
`examples/data/oswec.h5`, and a fresh output directory. It writes `flap.csv`,
`pto.csv`, and `wave.csv`; `source.npz` stores those three arrays under the
same names. `gate_spec.json` was written before inspecting this source run.
The [portable test](../../../test_oswec_regular_derived_source_ode4.py)
advances public `WEC.run(..., solver="wecsim_ode4")` using only HDF5 and
the derived wave settings, then checks all 601 output samples and four
PTO-work windows.

The largest measured source/Python pitch and pitch-speed gaps are
`1.05e-15` rad and `9.58e-16` rad/s. PTO torque differs by at most
`8.06e-9` N m. Full-run MATLAB and Python PTO work are both 5,028.9553 J;
their relative difference is `4.27e-14`. The predeclared limits are much
looser but distinguish the numerical route from the physical implicit
default: on this input the latter differs by `0.001138` rad pitch,
`13.14` N m torque, and `1.399%` PTO work. This demonstrates source-step
solver generalization without promoting the published physical-default
case or changing any badge.
