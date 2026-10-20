# RM3 PTO spring and pretension source trace

`source.npz` retains every 0.1 s sample of the 400 s [R2025b MATLAB
baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/38198224158)
for the pinned WEC-Sim RM3 example. The source checkout was
`0753b2e47f2457c078751dcfe5d251d1767b80ab`. Its PTO has 50 kN/m
stiffness, 5 kN pretension, and 1.2 MN s/m damping. The archive stores the
time grid; both bodies' surge, heave, and pitch positions and speeds; PTO
stroke, speed, force, and mechanical power; and the three exported PTO
parameters. It omits source motion from the Python solver inputs.

The fixture was made by selecting body CSV columns 1, 3, 5 and 7, 9, 11,
and PTO CSV columns 0, 3, 9, 15, 21 (zero-based), then calling
`numpy.savez_compressed`. The archived CSV SHA-256 values are:

| Source CSV | SHA-256 |
| --- | --- |
| `RM3_PTO_PRETENSION_Pretension_body1.csv` | `857dfb2d563836caa150ca3b93d3871de3dd196f83f0153d70898131741f7ed1` |
| `RM3_PTO_PRETENSION_Pretension_body2.csv` | `5bfb11a6e8f17000555d0e297e232e66c40d860eec823cc78cfe02850e02a14d` |
| `RM3_PTO_PRETENSION_Pretension_pto1.csv` | `ee965ad75d9d5fa151daee2858951abc5e2d027cb0c04cf4139a26fdda3b5238` |
| `RM3_PTO_PRETENSION_parameters.csv` | `731024e4df94828108b9fa2048be9b65b7a6831ff9c46fa703912a67799265bc` |

The bundled `examples/data/rm3.h5` SHA-256 is
`3bdfdcfa3dcd2db1851ab2026eb6b4fdc565ad86be8503e09b42231feee429fb`.
The test still accepts fresh MATLAB CSVs when they are supplied by the
reference-model workflow. This configured PTO is an RM3 example variant,
not an additional published input in the 63-case inventory.

The paired test also computes irreversible damper work as
`1_200_000 * pto_speed**2` from MATLAB's saved PTO speed. It gates separate
0–100 s, 100–250 s, 250–400 s, and full 400 s windows against the public
Python absorbed-power history at 1% relative error. The largest observed
window gap is 0.7934%.
