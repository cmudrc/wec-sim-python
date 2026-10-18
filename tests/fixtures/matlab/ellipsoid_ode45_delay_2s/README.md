# Published ellipsoid ode45 added-mass delay

`signals.npz` contains the first 2 s of the published
`Nonlinear_Hydro/ode45/Regular` R2025b source run. The pinned WEC-Sim core is
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications is
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The copied source model's
SHA-256 is `a5ebfaffdeaeb59b68665079ed1880c7657371f765db71e6c689f7cd2c4b989f`.
Its BEMIO-generated `ellipsoid.h5` SHA-256 is
`2eebf66e1ce1b22e229d59be05b2ebae648e675209c4c2469fc92eac5a7fc5ba`.

The local MATLAB R2025b Update 7 signal run changed only model `StopTime`
from 150 to 2 s. Simulink signal logging recorded the six acceleration
channels entering and leaving the active `10e-8` s `Transport Delay` under
`Rigid Body/Hydrodynamic Body`. At all 41 output samples, its input heave
acceleration equals the independently saved 150 s published body acceleration
exactly. The archived body and PTO force columns come from that 150 s run.
`tests/matlab/ellipsoid_ode45_delay_prefix.m` reproduces the signal log.

At each sample from 0.1 s, the delayed heave acceleration equals linear
extrapolation of the two preceding 0.05 s acceleration samples to
`t - 1e-7` s within `1.2e-17` m/s². The modified added-mass coefficient is
`-409527.37743276125` kg; multiplying it by the delay output reproduces the
source's logged applied added-mass force within `5.5e-12` N. Its adjusted
rigid mass is `742787.3346711756` kg, and the recorded heave force balance
closes within `6.9e-12` N.

These are source numerical identities, not an independent Python motion gate.
The published 150 s case remains qualified.
