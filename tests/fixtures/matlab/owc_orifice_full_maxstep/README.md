# OWC OrificeModel full-duration solver-step audit

Two independent local MATLAB R2025b Update 7 processes on Mali ran the
published 130 s `OWC/OrificeModel` input with pinned WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. Both use the same cleaned
seven-DOF HDF5 as the [published fixture](../owc_orifice_130s/README.md),
PM seed 128, `simu.dt=0.005` s, `ode23t`, and 26,001 output samples. The
scratch copies disable Simscape Explorer and post-run plotting and remove
unused `tf` and LPF50 declarations because the local installation lacks
Control System Toolbox. These edits are identical in both arms. Their
`wecSimInputFile.m` SHA-256 is
`90fb4b3c94ba524f1fefccb8044c1de4028f1183af65318337927e6da4dcc499`;
the unchanged `OWC_GBM.slx` and cleaned HDF5 hashes are
`53cdaf4830d7ac9340a17a8206b525b6c2b405ac18d8cc07912541ea2edec3cc`
and `fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`.

The nominal arm uses initialized `MaxStep=0.005` s. In the second arm, after
`initializeWecSim` and before `sim`, only the model's `MaxStep` changes to
`0.0025` s; `sim(..., simset('SrcWorkspace','parent'))` and `stopWecSim`
run in the same nested script context. Each source simulation took about
245 s. Both exported body, flexible mode, wave, component, and orifice
histories. The 500 components and all 26,001 wave samples match **exactly**;
both complete runs reproduce the corresponding earlier 10 s local audit at
every prefix body sample. The nominal and half-step `body.csv` SHA-256 values
are `acd8c28a155078da26d9abd3ef30fa1ae78c9e3b59aab5414d84d535c0724c`
and `2f60dabea3f56e331e589a0003d2090c18efbef4b06c7e704fe9903c297330c8`.
The shared wave and component CSV hashes are
`15b0f0e5ae2885ca0993831542a484c17ce93b63162e8542121a80474576216b`
and `a9effa6020e397f8913841de8572a7ff68810a9a82cb6419c6fa487e1df5d054`.
Against the archived hosted 130 s MATLAB run, the local nominal body differs
by at most `9.17` mm surge, `9.53` mm heave, and `0.252` mrad integrated
pitch despite identical waves and components. The step comparison therefore
uses the two same-platform Mali traces.

`source.npz` retains the common time and each arm's wave, components,
surge, heave, integrated unwrapped pitch, flexible displacement, and
output-grid orifice power in SI watts. Its SHA-256 is
`3d5bd33964de84e34812f5aaf01eac93cda4713e104eb013c06ca8128137873f`.
The pitch is the trapezoidal integral of saved pitch speed, as in the
published OWC fixture; raw Euler angles can jump branches near a singularity.

| Quantity | MATLAB MaxStep self-change | Independent Python gap vs nominal | Python gap vs half-step |
| --- | ---: | ---: | ---: |
| Surge peak / RMS | 34.84 / 14.84 mm | 166.01 / 131.32 mm | 155.67 / 118.61 mm |
| Heave peak / RMS | 8.63 / 2.81 mm | 39.89 / 20.01 mm | 37.36 / 19.97 mm |
| Unwrapped pitch peak / RMS | 5.03 / 0.855 mrad | 10.96 / 1.711 mrad | 5.93 / 0.897 mrad |
| Flexible displacement peak / RMS | 1.44 / 0.503 mm | 5.80 / 2.682 mm | 5.26 / 2.685 mm |

The independently advanced Python trajectory uses only the pinned HDF5 and
saved PM phases with public `WEC.run(..., orifice_force_path="published_owc")`.
Its wave elevation differs by at most `2.73e-14` m. Its PTO work is
`341,698.74` J versus `325,562.37` J nominal and `326,258.59` J half-step:
the MATLAB step change is `0.214%`, while the Python/half-step difference
remains `4.732%`. These measured residuals do not promote this derived
source arm to physical parity. MATLAB's full-duration
step change explains much of the pitch gap but only 18.4% of the squared
surge RMS gap; the physical coupled force-path issue remains separate. The
published OWC case stays qualified and badges are unchanged.

Refining only the independent Python integration from `0.005` to `0.0025` s
changes shared-sample surge, heave, pitch, and flexible displacement by at
most `0.0455` mm, `0.830` mm, `0.00963` mrad, and `0.116` mm. Python PTO work
changes by `0.0400%`. Against the same half-step MATLAB trace, the refined
Python surge RMS gap remains `118.59` mm, and its PTO-work gap is `4.774%`.
The ordinary Python test gates the Python step change without requiring the
remaining MATLAB/Python gap to stay large; a future model improvement can
reduce that gap without breaking this audit.

## Full-duration quarter MaxStep arm

A third local R2025b Update 7 process held `simu.dt=0.005` s, the published
130 s stop time, PM seed 128, HDF5, model, and all configuration fixed while
setting only initialized `ode23t` MaxStep to `0.00125` s. It completed in
226 s of wall time. The input-file SHA-256 is
`90fb4b3c94ba524f1fefccb8044c1de4028f1183af65318337927e6da4dcc499`;
the exported body and internal-orifice CSV hashes are
`9f4460a1b5ad579633ac987c8ad695924ef8fef067f95fd7cb811986718bd441`
and `8b0fb50d961038c8941e20769d4ef13b79faf69402fbc84fa1134b726ab2a22b`.
The compact `quarter_source.npz` SHA-256 is
`3a9c5d2cf8a104bd37c7a90553c2e405417af742afd34aad0d27d793e06df37a`.
Its first 30 s of body and flexible motion are byte-identical to the earlier
[MaxStep-only arm](../owc_orifice_dt_30s/README.md). All 500 PM components
and 26,001 incident-wave values equal the two previous local source arms
exactly. The orifice work integrates 427,235 strictly increasing internal
samples rather than the coarser output-grid record.

Relative to the 2.5 ms MaxStep source, another halving changes 130 s surge by
`133.15` mm peak and `82.24` mm RMS, and PTO work by `0.383%`. Against this
new source, the independently advanced public Python trajectory using the
opt-in `published_owc` force route differs by `77.19` mm peak and `49.55` mm
RMS surge, and `4.332%` PTO work. These improve on its comparison with the
2.5 ms arm but leave a material gap, while MATLAB's own surge self-change
has grown. The 130 s source has not established step convergence. This arm
remains a diagnostic; the published case and strict badge counts are
unchanged.
