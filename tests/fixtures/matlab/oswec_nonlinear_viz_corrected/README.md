# Corrected OSWEC nonlinear visualization baseline

This is the adopted `wec-sim-python` validation baseline for published
`Paraview_Visualization/OSWEC_NonLinear_Viz`. The pinned MATLAB source is
[WEC-Sim `0753b2e`](https://github.com/WEC-Sim/WEC-Sim/tree/0753b2e47f2457c078751dcfe5d251d1767b80ab)
with [Applications `d53d4d4`](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/d53d4d4c9eda2581f04204f5d394a6ef84bb099e),
run under R2025b. The [source workflow](../../../../.github/workflows/oswec-nonlinear-viz-corrected.yml)
and [MATLAB runner](../../../matlab/oswec_nonlinear_viz_corrected_baseline.m)
make disposable copies of the pinned case. The original source files remain
unchanged. Only the copied input's `simu.dt` changes from 0.1 to 0.0125 s,
refining the `ode4`, radiation, and nonlinear-force clocks together;
`simu.dtOut=0.1` s preserves the published 1,201-frame output and ParaView
grid. The 120 s duration, 40 s ramp, 2.5 m / 8 s regular wave, body and PTO
properties, nonlinear hydrodynamics, pressure logging, and ParaView setting
are held. A second copied arm uses `simu.dt=0.00625` s for a motion stability
check; pressure and scene logging are disabled in that arm.

The [hosted R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/39890460631)
produced the full corrected MAT histories and 3,603 wave/flap/base VTP files.
`tools/archive_oswec_nonlinear_viz_corrected.py` packed them without
rounding into `history.npz` (all 1,201 samples of both bodies' pose and
velocity, hinge PTO, and all three 1,042-facet pressure fields) and
`vtp_samples.zip` (the complete PVD and three source scene frames). The
source histories are checked again against freshly generated MATLAB output
in the hosted workflow. SHA-256:

| File | SHA-256 |
| --- | --- |
| `baseline.mat` from the run artifact | `0bb4d4381796060e1c67e9c5493fa6d2a3e0df55f682af18c2c807c3d634fdc8` |
| `stability.mat` from the run artifact | `3e939d4f1700c92cc0cc586ed8dc9a07f23aa67cbc538116a176c09fece94e70` |
| `history.npz` | `32e00658d9f0452b1cc2685d5f962c6304d541ead288526ab33446e9aa3080a1` |
| `vtp_samples.zip` | `e65a8f10ea7fb71574ef20cc93f3238ccf24b71be8a809d0d9cbc52896b0887d` |

The [ordinary Python gate](../../../test_oswec_nonlinear_viz_corrected.py)
advances the full 120 s physical flap motion from the bundled OSWEC HDF5,
STL, wave, hinge, and mass inputs. It does not prescribe MATLAB poses,
forces, or solver times. Against the 0.0125 s source, maximum active pose
gaps are 0.622 mm surge, 0.262 mm heave, and `0.0001351` rad pitch;
velocity gaps are 0.522 mm/s surge, 0.240 mm/s heave, and `0.0001095`
rad/s pitch speed. The 0.0125-to-0.00625 s source stability changes are
0.482 mm surge, 0.181 mm heave, and `0.0001014` rad pitch. The base remains
fixed. Hinge stroke and speed gaps match pitch and pitch-speed gaps; source
zero-PTO numerical force and power are below `2.6e-7` N m and `5.0e-8` W.

Facet pressures are recomputed from those independently advanced Python
poses. The maximum linear pressure gap is 1.10 Pa and the maximum physical
hydrostatic-plus-nonlinear pressure gap is 7.45 Pa. Eighteen of the
1,251,442 frame-facet positions switch wet/dry classification between the
two trajectories; each is within 1 mm of the free surface, where the
individual hydrostatic and nonlinear pressure terms jump in opposite
directions. The test checks this condition explicitly and gates their sum.
The complete 1,201-frame Python VTP scene is generated from Python motion
and pressure. Ordinary CI compares three archived source frames and all
3,603 PVD entries; the MATLAB workflow compares all 3,603 fresh source VTP
files. The largest flap vertex displacement is approximately 1.50 mm.

The unchanged published `simu.dt=0.1` s source trajectory differs from the
physical Python run by `0.007301` rad pitch. Its own change to the refined
0.0125 s source is `0.007192` rad. This correction establishes parity for
the documented repository-owned baseline; it does not claim physical parity
for the unchanged coarse-step source trajectory or require an upstream edit.
