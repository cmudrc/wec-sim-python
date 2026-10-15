# Published OSWEC hydraulic PTO layouts, 400 s

`source_50ms.npz` selects unrounded float64 MATLAB output for both published
`OSWEC_Hydraulic_PTO` and `OSWEC_Hydraulic_Crank_PTO` layouts from the
[passing R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37831809681).
It retains all 500 PM component entries and every fifth 10 ms output sample
through 400 s for wave elevation, flap motion and six-component excitation,
crank motion and torque, hydraulic pressures and flows, motor, and generator.
Each saved signed crank-work scalar is integrated from all 40,001 source
samples before selection. SHA-256:
`58c80c44eae64683da6d461b2942dd7c5fcf398fa06e48d9a49d4a389876b392`.

The test reuses `../oswec_multi_wave/oswec.h5` (SHA-256
`b9ccf3a32e0893c23b1645c15a2baef4a0c6f678b90d177d5797e3c65de4b7d3`).
Its frequency-domain coefficients are identical to the source run's BEMIO
HDF5. Four time-domain impulse-response arrays differ by at most `4.1e-11`
in stored units. The Python simulation advances every 10 ms step from the
pinned wave phases and HDF5, without reading any source motion, pressure,
controller, or force history as an input.

The hosted MATLAB source gate still checks every 10 ms output sample. This
portable gate checks every 50 ms output sample and the full-resolution work
integral, so it can run on each ordinary Python CI build.
