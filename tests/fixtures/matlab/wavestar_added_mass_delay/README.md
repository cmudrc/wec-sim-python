# WaveStar added-mass feedback audit

`source.npz` records a local MATLAB R2025b/Simulink run of the first 9.95 s
of the pinned 50 ms WaveStar NMPC application, before its controller activates.
It uses WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`,
WEC-Sim_Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`,
and `../wavestar_fixed_input_step/wavestar.h5`. A second run changes only
`body(1).adjMassFactor` from its source default of 2 to 1 in a disposable
application copy. Both runs shorten the end time and disable post-run plotting;
the wave, geometry, solver step, and controller settings are unchanged.

The default prefix reproduces the archived published source pitch and speed
within floating-point roundoff. Its float rigid mass is 3.075 kg, while the
source passes 10.7767874 kg to Simscape after shifting twice the trace of the
translational infinite-frequency added mass into the body. The compensating
added-mass matrix is applied through the source library's 100 ns transport
delay. The saved force starts at zero; at 50 ms it equals that matrix times
the initial acceleration. At later logged times it closely follows a linear
extrapolation of the two preceding acceleration samples. With factor 1, the
initial pitch acceleration changes from 0.338370 to 0.519208 rad/s² and the
first-9.95-s pitch history changes by 0.00272635 rad.

The default run also saves the BC rod and motor's six-component positions,
velocities, and accelerations at all 200 samples. Their body origins start at
the world origin. The rod rotates with the B-to-C direction and translates so
its C joint follows the arm; the motor origin rotates about B. A separate
portable kinematics gate compares this geometry with the saved source poses
and velocities before those bodies' mass and gravity enter the Python solve.

This is a source numerical diagnostic, not proof that the published trajectory
is physically converged or that the Python model should adopt the delayed
feedback. It does not close the 225 s WaveStar parity case.

To regenerate, run `tests/matlab/wavestar_added_mass_delay_probe.m` twice with
the local pinned core, applications, and this Python checkout as its first
three arguments. Use fresh output directories; pass `1` as a fifth argument
for the alternate run. Each run writes `added_mass_probe.mat`. Convert those
files with SciPy:

```python
from pathlib import Path
import numpy as np
from scipy.io import loadmat

default = loadmat(Path('/tmp/default') / 'added_mass_probe.mat')
factor_one = loadmat(Path('/tmp/factor1') / 'added_mass_probe.mat')
out = Path('tests/fixtures/matlab/wavestar_added_mass_delay/source.npz')
np.savez_compressed(
    out,
    time=default['time'].ravel(),
    pitch=default['pitch'].ravel(),
    speed=default['speed'].ravel(),
    acceleration=default['acceleration'],
    applied_added_mass=default['appliedAddedMass'],
    modified_added_mass=default['modifiedAddedMass'],
    original_added_mass=default['originalAddedMass'],
    adjusted_mass=default['adjustedMass'].item(),
    rigid_mass=default['rigidMass'].item(),
    mass_factor=default['massFactor'].item(),
    factor_one_pitch=factor_one['pitch'].ravel(),
    factor_one_speed=factor_one['speed'].ravel(),
    factor_one_acceleration=factor_one['acceleration'],
    aux4_position=default['aux4Position'],
    aux4_velocity=default['aux4Velocity'],
    aux4_acceleration=default['aux4Acceleration'],
    aux5_position=default['aux5Position'],
    aux5_velocity=default['aux5Velocity'],
    aux5_acceleration=default['aux5Acceleration'],
)
```

The portable Python test checks source input identity, the mass shift, the
applied-force timing, and the alternate-factor sensitivity.
