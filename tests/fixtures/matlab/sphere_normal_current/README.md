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

The test advances both Python trajectories independently from the wave and
current inputs, and retains the original position, speed, current-effect, and
source-force limits. It also checks that the source control and current waves
are identical. Saved MATLAB states are used only by the separate force-law
diagnostic, not to drive the coupled Python motion. This is a derived case and
does not change a published-case badge.
