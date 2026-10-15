# Local OWC OrificeModel source prefix

`source.npz` records a bounded source-to-source check, not a Python parity
baseline. The archived published 130 s source output came from [GitHub Actions
run 37799665808](https://github.com/cmudrc/wec-sim-python/actions/runs/37799665808)
at commit `5dbf0d5fc228bcbd39275ec8c267ef0e43289bb2` on an Ubuntu runner.
The local run used MATLAB R2025b Update 7, WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab`, and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Both used the exact BEMIO
HDF5 in `../owc_orifice.h5` (SHA-256
`fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`).
The fixture SHA-256 is
`9dc359ff2ef46dc298f3a753c3eb4cf08415d17bfefab1c275682ac7098e6086`.

The local source was a disposable copy of the pinned `OWC/OrificeModel`. The
6 s run changed `simu.endTime` from 130 to 6, disabled Simscape Explorer,
removed the unused `tf`/LPF50 declarations because this MATLAB installation
lacks Control System Toolbox, and replaced post-run plotting with a no-op.
Inspection of `OWC_GBM.slx` found no reference to the removed filter variables.
No physical model setting, solver, output step, wave seed, or orifice parameter
was changed. A second local 6 s run yielded byte-identical exported body,
wave, and orifice CSVs. Changing only `simu.endTime` from 6 to 10 in the same
scratch model yielded a byte-identical first 6 s of those outputs and the wave
components. The local 6 s result is thus repeatable, and truncating it to 6 s
does not explain its difference from the archived run. Other host or source
execution differences have not been isolated.

The `source.npz` arrays are `body` (time plus six rigid positions and six
velocities), `flex` (saved flexible-mode output), `wave` (time and elevation),
`components` (frequency, amplitude, bin width, phase), `orifice` (local
internal evaluations), and `archived_orifice` (archived evaluations through
6 s). The archived rigid/flexible histories and input data remain in
`../owc_orifice_prefix.npz`.

Before the local run, the gates for local versus archived output were fixed at
`1e-10` absolute for the incident wave and all 500 components and `1e-5`
absolute for body position, body speed, flexible displacement, and flexible
speed. The wave and components were bit-equal; rigid motion matched within
`1e-12` through 2 s. The four motion gates failed:

| Quantity | Maximum absolute gap | Time |
| --- | ---: | ---: |
| Rigid position (surge) | 0.001054317 m | 6 s |
| Rigid speed (surge) | 0.000701002 m/s | 4.81 s |
| Flexible displacement | 0.0000337022 m | 5.20 s |
| Flexible speed | 0.000332673 m/s | 4.505 s |

The archived run logged 1,639 internal orifice evaluations through 6 s;
the local run logged 1,437. At 6 s their reaction forces were -65.370925 N
and -64.405790 N, respectively. The fixture's ordinary Python test keeps
this failed source-to-source gate visible. It does not promote
`OWC/OrificeModel` to a published-input gate or prove which host/model factor
caused the difference. The physical coupled Python path and explicit published
Simulink force route retain their separate limits in `PARITY.md`.
