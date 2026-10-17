# Moving regular Morison source force

`MORISON_MOVING_REGULAR_states.csv` is the eight-state option-1 force export
from [R2025b run 37768083149](https://github.com/cmudrc/wec-sim-python/actions/runs/37768083149),
using WEC-Sim core commit `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
`tests/matlab/reference_model_baseline.m` generates the CSV in the
`MORISON_MOVING_REGULAR_states.csv` block. Each row contains time, six body
position components, six velocities, six accelerations, and six source force
and moment components. The archived CSV SHA-256 is
`466ff5a912ca07062d699f9275e6ea2db2ae42b20ee3b6844a967e53f96ca827`.

The portable test evaluates Python's source-law diagnostic at all eight saved
states, including a wave ramp, nonzero rotation, and an emerged element.
These prescribed states do not establish coupled moving-body parity. The
separate physical axial-law test protects its rigid rotation and dissipative
relative-velocity behavior; [WEC-Sim issue #1664](https://github.com/WEC-Sim/WEC-Sim/issues/1664)
tracks the nonorthogonal matrix retained only by the Python source diagnostic.
