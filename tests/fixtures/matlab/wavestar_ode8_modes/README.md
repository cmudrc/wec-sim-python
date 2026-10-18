# WaveStar radiation-mode solver probe

`amplification.npz` pairs the 92 active state-space radiation poles in the
pinned `../wavestar_fixed_input_step/wavestar.h5` with their one-step
amplification under MATLAB R2025b Simulink fixed-step `ode8` at `dt=0.05 s`.
The float's surge, heave, and pitch inputs select the active blocks. The
probe is a block-diagonal linear model with no forcing, controller, or body
coupling; it measures only the solver's treatment of those poles. It cannot
establish closed-loop WaveStar motion parity.

To regenerate from this checkout, first save the poles for MATLAB:

```bash
python - <<'PY'
from runpy import run_path
from scipy.io import savemat
sweep = run_path('tests/test_wavestar_ode8_modal_stability.py')
savemat('/tmp/wavestar_modes.mat',
        {'eigenvalues': sweep['_active_radiation_poles']()})
PY
```

Run `tests/matlab/wavestar_ode8_modal_probe.m` in MATLAB R2025b with
`/tmp/wavestar_modes.mat` as `modesFile` and a temporary output path as
`outputFile`, then convert its arrays:

```bash
matlab -batch "addpath(fullfile(pwd,'tests','matlab')); wavestar_ode8_modal_probe('/tmp/wavestar_modes.mat','/tmp/wavestar_ode8.mat')"
python - <<'PY'
from scipy.io import loadmat
from numpy import savez_compressed
source = loadmat('/tmp/wavestar_ode8.mat')
savez_compressed('tests/fixtures/matlab/wavestar_ode8_modes/amplification.npz',
                 eigenvalues=source['eigenvalues'].ravel(),
                 ode8_amplification=source['amplification'].ravel())
PY
```

The portable Python test re-derives the poles from the HDF5 and compares
the archived source amplification with RK4, DOP853, and exact linear flow.
