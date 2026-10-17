# Sphere normal/tangential Morison current pair

`source.npz` preserves the full 4,001-sample no-current control and current-run
body histories, the direct MATLAB current-force increment, and both wave
histories from [R2025b run 39068743962](https://github.com/cmudrc/wec-sim-python/actions/runs/39068743962).
The run used pinned WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The current
run changes only the derived normal/tangential Sphere case by adding a
0.8 m/s one-seventh-power current. The two saved wave CSVs are byte-identical.

The unrounded NPZ SHA-256 is
`25d29ed94c69fc61010b6e1f8c8671ff4784c85064dd323a89d37c229065a4f7`.
The source body CSV SHA-256 values are
`767820ea5d681862a515dd1d20382fa8abb330c1cc20941cb5687a13ba244733`
without current and
`dca7a029ead08e94462347a996a26a8e7044792badf3860fa8ff336220fe4e23`
with current. The direct-force increment CSV SHA-256 is
`6ca19418233e2562bd73acef8130f3d86dedea31042583512513c60e781bf81f`.
The source HDF5 container SHA-256 is
`ad78ad6a105806d576aa27aaae95fe3c0d98ca94dd16bf52fa58ae9d20abc81d`;
all 46 datasets, including NaN locations, are exactly equal to the bundled
`examples/data/sphere.h5`. The portable test therefore uses that bundled file.

The portable normal/tangential and current tests advance the Python body
independently from the wave and current inputs. They retain the original
full-duration position, speed, Morison-force, and current-effect limits. The
current test also checks that the source control and current waves are
identical. Saved MATLAB states are used only by the separate force-law
diagnostics, not to drive coupled Python motion. These are derived cases and
do not change published-case badges.

## Local source reproduction

On 2026-10-17, R2025b on macOS ARM64 reran both source arms from disposable
copies of the pinned `Free_Decay/1m-ME` input and Sphere model, using the exact
HDF5 from the archived artifact. The input file was reset to its pinned state
between `reference_model_baseline('SPHERE_MOVING_MORISON_NORMAL')` and
`reference_model_baseline('SPHERE_MOVING_MORISON_NORMAL_CURRENT')`.
The no-current 4,001-row body CSV and both wave CSVs are byte-identical to
the hosted artifact. Across all 37 logged current-body columns, the largest
local-to-hosted numerical difference is `1.001e-11` in its logged units; the
seven-column direct current-force increment differs by at most `1.019e-10`.
Both paired Python tests also pass against the local source CSVs. This checks
the portability of the pinned MATLAB source data; it does not add a new case.
