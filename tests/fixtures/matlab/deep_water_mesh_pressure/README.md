# Deep-water regular-wave facet pressure reference

`source.csv` was generated locally with MATLAB R2025b Update 7 from pinned
WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab`, using
[`deep_water_mesh_pressure_reference.m`](../../../matlab/deep_water_mesh_pressure_reference.m)
and the source `nonFKForce` function. The CSV SHA-256 is
`25783a1647273b7023f4424ba896410a4306563830f563ed5587279d0276f28e`.

Two upward-facing 0.5 m² facets have body-local centers at `z=-2` and
`z=-0.25` m. The three rows prescribe a regular-wave crest at zero heave,
a trough at zero heave, and a crest after 1 m upward heave. The shallow
facet is dry in the last two rows. The 13 columns are time, heave, wave
elevation, two nonlinear pressures, two mean-surface linear pressures, and
six nonlinear-minus-linear Froude-Krylov wrench components.

The portable test checks pressures and the correction wrench to `1e-8`
Pa or N, and separately checks hydrostatic pressure against the analytic
facet depths. This is a prescribed-pose force-law gate for deep water. It
does not establish a published input's independent body-motion parity or
change a reference-model badge.
