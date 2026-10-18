# WaveStar AR forecast on independently advanced motion

`python_motion.npz` saves the 4,501 pitch and requested-torque samples from
the autonomous Python 225 s `run_wavestar_published(...,
solver="wecsim_ode8")` run used for the [independent QP input
gate](../wavestar_nmpc_independent_qp/README.md). It used the pinned
`wavestar_fixed_input_step/wavestar.h5` fit and the 500 published incident
wave phases; MATLAB body motion, controller states, and torque requests were
not runtime inputs. The portable test rebuilds the observer's excitation
estimate from these Python pitch and command samples and requires its active
states and requests to agree with that independent QP fixture.

`python_moment.csv` records those 4,501 reconstructed excitation estimates
at the 0.05 s controller clock. On 2026-10-18, MATLAB R2025b Update 7 on
mali ran [`wavestar_nmpc_independent_predictor_replay.m`](../../../matlab/wavestar_nmpc_independent_predictor_replay.m)
using the pinned WaveStar controller's 180-sample, order-18 `ar` fit and
40-sample `forecast` call. `matlab_forecast.csv` contains all 216 refits at
integer seconds from 10 through 225 s. MATLAB also called the pinned
`predictor.m` on every 0.05 s input from 10 through 225 s, preserving its
persistent fit between refits. `matlab_active_forecast.mat` contains every
forecast consumed by the controller from 15 through 225 s: 4,201 times 40
horizon values. The underlying Applications commit is
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`; its `predictor.m` has
SHA-256 `35cd4e524641de6281a9bd7d5ac3316832780db6dfc9213135aa27f577b903c2`.

SHA-256 of `python_motion.npz`, `python_moment.csv`,
`matlab_forecast.csv`, and `matlab_active_forecast.mat`, respectively:

```text
2abec603f3d0975d37263fb32622357f7a05edc1694a1fa2edb121df010dd5dc
56c46a6be44863ff29444bcc70d2920dce989aebca9df66d80ab52f2b3bfb0e0
374041ef04795b46cf98e7a05ae4111ac81a5e91ca163078d413e448e05fb530
95df40f7cc6ca4468f20cd5d01ff696040da9e3bc7d17f0641a85e854a0fd5f0
```

The Python predictor differs from direct MATLAB forecasts by at most
`5.77560e-8` N m over all 8,640 predicted values, at the first 10 s refit;
RMS difference is `2.06425e-9` N m. Later refits differ by at most
`3.49587e-12` N m. The portable limits are `1e-6` N m maximum and
`1e-7` N m RMS. Across all 168,040 active-controller forecasts, the
maximum and RMS differences are `4.55458e-12` and `3.85620e-13` N m,
below the same limits. A bounded precontrol probe found an off-refit
`0.09197` N m difference from 10 to 11 s, during the first fitted model's
hold interval; the predictor refits before the controller starts at 15 s.
This validates the predictor component on independent Python motion inputs.
It does not advance MATLAB and Python closed-loop plants from the same
initial state, so the published 225 s motion remains open.

To regenerate the MATLAB forecast CSV, first export the observer excitation
estimate from the autonomous Python run to `python_moment.csv` in the time
order checked by the portable test. Then run
`wavestar_nmpc_independent_predictor_replay(applicationsRoot)` from local
MATLAB with `tests/matlab` on its path, where `applicationsRoot` is the
pinned WEC-Sim_Applications checkout. The script reads the committed input
CSV and overwrites the two MATLAB output fixtures; rerun the portable test
afterward. The MAT-file header may change its byte hash when regenerated;
the test gates numeric values.
