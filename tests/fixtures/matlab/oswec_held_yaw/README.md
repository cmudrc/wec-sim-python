# Published OSWEC held-yaw source traces

`source.npz` retains float64 arrays from the successful R2025b
[PassiveYawRegression source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38403345483)
at pinned WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Seeds 1, 2, and 3 each
cover the published 250 s at 0.01 s. The archive preserves the source's
500-component frequency, spectral amplitude, interval, and phase arrays;
wave elevation; flap yaw and speed; yaw excitation; PTO speed, torque, and
power. Values are taken directly from the artifact CSV columns without
resampling or rounding. SHA-256: `129977c8d693d3e3d12ed4cee8ddd3b22e2746d7a3bd03d901a5d11348697095`.

The test reuses [`../oswec_multi_wave/oswec.h5`](../oswec_multi_wave/oswec.h5).
All 68 HDF5 datasets match the source artifact's generated `oswec.h5`
elementwise, including NaN positions. The files differ in HDF5 storage
representation. Python generates wave phases from the published MATLAB
substream seeds; archived phases are comparison outputs, not inputs.

These gates reproduce the published `ode4` numerical path through the opt-in
`wecsim_ode4` solver. The physical implicit-radiation default and source
step sensitivity remain separate, so this fixture does not add a published
physical-motion badge.

## Repository-owned continuous-heading control

[`continuous_heading.patch`](continuous_heading.patch) changes only the two
body yaw thresholds from 1° to zero in the pinned Applications input above.
The [three-seed source and physical gate](https://github.com/cmudrc/wec-sim-python/actions/workflows/oswec-passive-yaw-continuous-step.yml)
checks out that unmodified upstream revision and applies the patch in CI.
WEC-Sim's pinned `irregnLYaw.m` uses threshold zero to interpolate heading
coefficients every step. This control has no discrete held-heading events.
It is a derived input owned by this repository; it does not change the
published Applications input or promote the original 1° case.
