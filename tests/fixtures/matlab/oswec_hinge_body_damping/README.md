# OSWEC hinged-flap body damping

`source.npz` contains local MATLAB R2025b Update 7 output from the pinned
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
It starts from Applications `Nonhydro_Body/wecSimInputFile.m` and the core
`examples/OSWEC/OSWEC.slx`. Both arms shorten the run from 400 to 150 s.
The damped arm additionally sets only
`body(1).linearDamping = diag([100000 0 0 0 2000000 0])`;
the PTO remains undamped. The two solver-only arms retain all WEC-Sim
0.1 s clocks and change only the Simulink `ode4` fixed step to 0.025 s.
`tests/matlab/oswec_hinge_body_damping_reference.m` exports the trajectories
from disposable application copies using the pinned OSWEC HDF5 in
`tests/fixtures/matlab/oswec_multi_wave/oswec.h5`.

The unchanged local 0.1 s arm reproduces the archived published flap and
base traces within `1e-15`. MATLAB's body-damping output equals `D @ v`
on all six saved flap velocity components within `1.1e-9` N or N m; that
logged term is subtracted in the force balance. The damper changes source
flap pitch by up to `0.03365` rad over 150 s. The public Python implicit
solver independently advances both arms from the HDF5 and regular wave.
Against coarse MATLAB, maximum damped flap errors are `18.71` mm surge,
`6.85` mm heave, `0.00398` rad pitch, and `0.00301` rad/s pitch speed;
the published OSWEC absolute motion limits pass.

The tighter predeclared change-from-control gate **fails** against 0.1 s
MATLAB: surge and pitch changes differ by `3.83` mm and `0.000830` rad,
above `1` mm and `0.00015` rad limits. With only MATLAB's solver step
refined to 0.025 s, those differences fall to `0.230` mm and
`0.0000500` rad; all six unchanged differential position and speed gates
pass. The opt-in public `wecsim_ode4` path pairs both coarse trajectories
to below `5e-14` rad pitch. This supports the body damping force law and
separates a coarse solver-step effect. It is a derived configuration of an
already counted published case; the case badges do not change.

SHA-256 of `source.npz`:
`7a1cd4776883c050572d3b0bcd81f633934c52d694e2905b887f091e6a3b9072`.
The disposable control and damped input files have SHA-256 values
`9e1d0ce778ce2970eeb03518fb5da941dd00179a1d46d1456db1ab17d278d9fa`
and `7fbf897a493026bc4c15c0f5455520063899b8bea1e73e2c63ca46c29f328d1f`.
