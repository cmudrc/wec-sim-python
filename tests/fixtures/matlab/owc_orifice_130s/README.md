# Published OWC OrificeModel source trace

`source.npz` is a compact copy of outputs from the successful [130 s R2025b
baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37799665808),
artifact `matlab-OWC_ORIFICE-baseline`, run commit
`5dbf0d5fc228bcbd39275ec8c267ef0e43289bb2`. The workflow pins WEC-Sim
core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The archived
`test17a_clean.h5` is byte-identical to `../owc_orifice.h5` (SHA-256
`fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`).

The fixture preserves all 26,001 output times, wave elevations, selected rigid
positions and speeds, flexible displacement and speed, all 500 source phases,
and the five orifice parameters as float64. `body_position[:, 2]` is the
trapezoidal integral of saved pitch speed, matching the existing fresh-source
gate's unwrapped-pitch comparison. Orifice force and numeric kW power from the
418,828 internal evaluations are interpolated onto the output grid; power is
converted to SI watts. The source files have these SHA-256 hashes:

| Archived source file suffix | SHA-256 |
| --- | --- |
| `components.csv` | `a9effa6020e397f8913841de8572a7ff68810a9a82cb6419c6fa487e1df5d054` |
| `OrificeModel_body1.csv` | `8153820cc0de55b347ca0f7330eea49d1e5c9e71176c8145898c5687b97219b2` |
| `orifice.csv` | `22314e78887134feebc28f2a54e1e5157d39769ce534229367f44dc22cd43096` |
| `Flex_out.mat` | `825ff26839df1fccb5c48533ec71511671240d5693450cf0cb99df1518c6e20c` |
| `wave.csv` | `15b0f0e5ae2885ca0993831542a484c17ce93b63162e8542121a80474576216b` |
| `parameters.csv` | `ab0b29d67f2651a3bf48a41889362b1d15790ad569e2ab96aff218a9b92d43db` |

The packed fixture's SHA-256 is
`149af08ecdc87a13c072149957aa8b60f577fb639eae676d8e5b2dfa55ff2630`.
`tests/test_owc_orifice_130s_portable.py` uses only the saved input phases and
HDF5 to independently advance the public `WEC.run` path over 130 s. It keeps
the existing fresh-source motion, PTO, work, and first-ten-second limits and
also checks incident-wave elevation. The opt-in `published_owc` force route
follows the source Simulink wiring. Its passing comparison is a
**source-convention gate**; the physical coupled default still diverges late,
and source operation above the orifice Mach threshold is not validated as
physically accurate. The case remains `qualified_motion` and gains no badge.
