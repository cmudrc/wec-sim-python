# Unchanged RM3 PTO Extension physical comparison

This fixture archives the **unchanged published** 30 s, 0.1 s `ode4`
`RM3_PTO_Extension/float` and `/spar` source runs. The six arrays in
`source.npz` are the complete 301-row, 25-column MATLAB body and PTO CSVs,
loaded as `float64` without interpolation or rounding. The exact generated
`rm3.h5` used by those runs is also archived. Both inputs have no wave, zero
PTO spring and damping, and a 5 m initial PTO extension.

- Source: [successful R2025b MATLAB reference run 38282150325](https://github.com/cmudrc/wec-sim-python/actions/runs/38282150325), repository commit `7c403d3e2e07b9940e4a2a280d375000b8985fea`.
- Pinned WEC-Sim core: `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
- Pinned Applications: `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
- `rm3.h5` SHA-256: `dac35aa42321e1ade04e8935b38680d7b70e7b247edac9b5b800d01142fcc259`.
- `source.npz` SHA-256: `535eaabbf0ec68949ffffff215548003696a0ad99029214a396b37d6fa67e00c`.

The nine time, active heave/speed, and PTO columns match the repository's
earlier `rm3_pto_extension_step.npz` published 0.1 s arms exactly. The source
HDF5's active center-of-gravity, volume, restoring, added-mass, and radiation
IRF datasets equal those in `examples/data/rm3_end_stops.h5` by value. The
portable test checks both identities and can compare all six arrays to a
fresh MATLAB reference output when `WEC_SIM_MATLAB_MODEL_OUTPUT_DIR` is set.
The [dedicated source workflow](../../../../.github/workflows/rm3-pto-extension-published-physical.yml)
regenerates both unchanged inputs from the pinned checkouts and runs that
fresh comparison together with the independent Python gate.

`tests/test_rm3_pto_extension_published_physical.py` advances a public
`run_case` using physical implicit added mass, the exact source HDF5, no
wave, the original 0.1 s grid and 30 s duration, and the published initial
coordinates. It then checks all six positions and velocities of both bodies,
PTO stroke and speed, and zero force and power against the complete MATLAB
records. The limits are **inherited unchanged** from
`tests/test_rm3_pto_extension_parity.py`: active float heave/stroke 90 mm
and speed 210 mm/s; active spar 18 mm and 8 mm/s; the original limits on
all other channels also apply. The physical run receives no MATLAB states,
accelerations, forces, or solver steps. An opt-in source-`ode4` run is used
only as a negative control: it nearly reproduces the published numerical
history, whereas the physical trajectory is measurably distinct.

On local arm64, the physical float's largest heave/stroke gap is 85.849 mm
and speed gap 202.600 mm/s, leaving margins of 4.151 mm and 7.400 mm/s.
The spar gaps are 14.877 mm and 6.522 mm/s, leaving margins of 3.123 mm
and 1.478 mm/s. These are formal passes under the inherited gates with
limited margin. The original limits preceded this source-step investigation
but were not prospective thresholds chosen blind to all earlier errors.

The published source remains numerically sensitive. Changing only its
Simulink `FixedStep` from 0.1 to 0.0125 s removes over 99% of the float's
and over 99.99% of the spar's squared heave/speed gap to a refined source
run; see `../rm3_pto_extension_fixedstep_only/README.md`. This test therefore
does not claim that the coarse published trajectory is source-converged or
that Python reproduces it pointwise. The separate corrected fine-step
baseline and the original published comparison retain distinct fixtures.
Case status and badge interpretation require separate review.
