# OWC OrificeModel 10 s internal-step audit

`source.npz` stores two local MATLAB R2025b Update 7 source runs on all 2,001
0.005 s output samples. Both use pinned WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab`, Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and the seven-DOF HDF5
`../owc_orifice.h5` (SHA-256
`fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`).
The compressed fixture SHA-256 is
`d140edad47e60f230a72d7ba3cfe3628aa604bdddcdb99671b4aeb5587d03d35`.

The scratch application changes only `simu.endTime=10`, disables Simscape
Explorer, removes unused `tf`/LPF50 declarations because this installation
lacks Control System Toolbox, and disables post-run plotting. These same
changes are in both runs. The first uses the initialized `ode23t`
`MaxStep=0.005` s. The second changes only model `MaxStep` to `0.0025` s
after `initializeWecSim`, then invokes `sim(..., simset('SrcWorkspace',
'parent'))` and `stopWecSim` from the same nested script context as `wecSim`.
The output grid, `simu.dt`, wave seed 128, HDF5, solver, and orifice parameters
stay fixed. The source model file itself is not saved with a changed MaxStep.

The arrays `body_nominal` and `body_half` contain time, six rigid positions,
and six velocities. `flex_*` contain the saved ten-channel `Flex_out` values;
columns 0 and 1 are flexible displacement and speed. `wave_*` and
`components_*` store incident elevation and all 500 PM components. The wave
and component arrays are bit-identical across the two runs. The original
source exports had SHA-256 hashes `79ac8317c4c7feca8e61a748a7d7cfc2d7f3d515a574757e645ef1ed6f378608`
and `a3ffba15c25126f2cc3a22d84d3f259f814e6a3b95b26d4834251ced84e61593`
for the nominal and half-step body CSVs.

Before the run, motion-stability limits were set at `1e-4` absolute for body
position/orientation, body velocity, flexible displacement, and flexible
speed; wave and component limits were `1e-10` absolute. The wave and
component limits pass exactly, but every active motion limit fails:

| Quantity | Maximum MATLAB self-change |
| --- | ---: |
| Rigid surge | 6.141 mm |
| Rigid heave | 0.586 mm |
| Relative 3D rotation angle | 0.005850 rad |
| Rigid surge speed | 2.526 mm/s |
| Flexible displacement | 0.1915 mm |
| Flexible speed | 2.047 mm/s |

Raw Euler roll/yaw sometimes changes branch by π near the pitch singularity;
the rotation row is the physical relative angle from the saved XYZ Euler
triplets. This local 10 s result is not a convergence proof or a hosted
published 130 s source run.

The ordinary Python test also independently advances the public `WEC.run`
model from the pinned phase seed and HDF5, with no source motion or force as
input. Against the half-step source, the explicit `published_owc` force route
has maximum surge/heave/pitch/flexible-position gaps of
6.691 mm / 1.063 mm / 0.005929 rad / 0.2421 mm. Its regression screens of
8 mm / 1.5 mm / 0.007 rad / 0.3 mm were set after this exploratory run; they
are diagnostic, not published-case parity gates. The physical coupled path
still differs by 82.15 mm heave and 20.94 mm flexible displacement. Thus
solver-step sensitivity helps explain surge and pitch differences, while it
does not remove the larger flexible-force-routing gap. See
[`PARITY.md`](../../../../PARITY.md).
