# Local WaveStar NMPC source diagnostics

The pinned WaveStar application uses System Identification Toolbox `ar` and
`forecast`, which are unavailable in the MATLAB installation on mali. The two
functions in this directory are a **test-only** replacement for the
application's numeric, univariate AR(18) predictor. `ar.m` fits the same
forward-backward modified-covariance regression used by the published case;
`forecast.m` advances its 40-step recursion. They are not general replacements
for System Identification Toolbox objects or APIs. Put this directory on the
MATLAB path only for a disposable application copy.

For a local source run, use pinned WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Copy the Applications
`WECCCOMP` directory outside both checkouts. Place
[`wavestar.h5`](../../fixtures/matlab/wavestar_fixed_input_step/wavestar.h5)
in that copy's `hydroData` directory. Disable only its
`WECCCOMP_Nonlinear_Model_Predictive/userDefinedFunctions.m` plotting script;
it expects an interactive workspace variable and is called after simulation.
From the copied `WECCCOMP_Nonlinear_Model_Predictive` directory, add the pinned
core's `source` tree and this directory to the MATLAB path, then call `wecSim`
as a script. The application keeps its published controller, sea, `ode8`
solver, and 225 s
duration. WEC-Sim writes `output/WaveStar_matlabWorkspace.mat`.

This route is a **diagnostic source substitution**, not an unchanged published
MATLAB run. The archived unchanged-source baseline and its independent Python
checks remain authoritative. The temporary predictor must be compared to that
baseline before a local trajectory is used to investigate a new mechanism.
Download the `matlab-WECCCOMP_NMPC_SOURCE-baseline` artifact from the linked
source run. With its `WECCCOMP_NMPC_SOURCE_controller.mat` and the pinned
application folder, call `check_archived_forecast(controllerFile,
applicationFolder)` to verify the replacement on all 4,501 estimator states.

Against the archived [published source artifact](https://github.com/cmudrc/wec-sim-python/actions/runs/37876181993),
replaying the helper on all 4,501 saved estimator states changes forecasts
after 11 s by at most `2.02e-11` N m. The unstable first fit differs by
`0.0525` N m at most against a `5.24e6` N m source forecast. A disposable
225 s run with the helper reproduces the 500 wave components and all 4,501
wave samples exactly at stored precision. Its body and PTO columns match the
archived first 20 s within `1.5e-7` in their logged units. Feedback amplifies
small later differences: the full run reaches `0.01412` rad pitch and
`0.882` N m command differences from the unchanged hosted source. Those are
**source-to-source differences**, not Python errors or a parity gate.
