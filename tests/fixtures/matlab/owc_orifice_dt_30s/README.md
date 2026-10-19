# OWC OrificeModel 30 s source-clock audit

Three local MATLAB R2025b Update 7 runs on Mali compare the pinned
`OWC/OrificeModel` with `simu.dt` and `ode23t` `MaxStep` both set to 5, 2.5,
and 1.25 ms. The first arm is the first 30 s of the archived 130 s nominal
run; the two shorter runs change only `simu.dt` and the stop time in disposable
copies. All use PM seed 128, the same cleaned seven-DOF HDF5, 10 s ramp,
15 s radiation memory, and the published geometry, damping, and orifice.
Simscape Explorer and unused Control System Toolbox declarations were removed
identically in all arms, as described in the [full-duration audit](../owc_orifice_full_maxstep/README.md).

Core WEC-Sim is `0753b2e47f2457c078751dcfe5d251d1767b80ab`;
Applications is `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The
unchanged model and HDF5 SHA-256 values are
`53cdaf4830d7ac9340a17a8206b525b6c2b405ac18d8cc07912541ea2edec3cc`
and `fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`.
The 2.5 and 1.25 ms input-file hashes are
`eeb18fb92cf41da204c39f57f5ff871d1e7ca3542e90873f3f4c3a3e9cceacae`
and `885146629fa94b7093068a078049d798c628b5a5aa52468fd4f71529d337120a`.
Their full body CSV hashes are
`16ad049b6d7063736fd95593cdc136964c355a078c03af4492f366b9d7148b02`
and `5a34481d371869802c6d575329107e07b2014c88a9fb627a514d61d40d485ba9`.
The compact `source.npz` SHA-256 is
`6b7418b5fe8e508c0115f527c9ccfb2d3932950122f38f4f6ec8ede800a1db51`.
The separate `internal_trace.npz` SHA-256 is
`c949f53ea2b75be5d278794f3454e857e89fcfe582c8114688337ae6f4c74453`.
One additional 30 s arm keeps `simu.dt=0.005` s and changes only initialized
`ode23t` MaxStep to `0.00125` s. Its input-file, body CSV, and compact
`maxstep_only.npz` SHA-256 values are
`8417681e1036c2d4f67e87fa4dc35c47a819022079b98012020bf015ebb8eb18`,
`b6f4171c4b2ceb39f42977618689ef150408bef23b8eee603d81f113cc28912a`,
and `b2f6aac77e08bb9e0fa2f5721cf356ef9e964396cbad14ce9d0e55ba5f360370`.

The fixture retains all three source trajectories on their common 5 ms grid,
the exact 500 PM components and incident wave, and PTO work integrated on
each arm's strictly increasing internal orifice timestamps. Its motion
columns are surge, heave, unwrapped pitch, and flexible displacement;
velocity columns follow the same order. Unwrapped pitch integrates saved
pitch speed, avoiding Euler-angle branch jumps.

| Source step | 30 s surge shift from preceding arm | Heave shift | Unwrapped pitch shift | Internal PTO work |
| --- | ---: | ---: | ---: | ---: |
| 5 ms | — | — | — | 56.184 kJ |
| 2.5 ms | 10.151 mm | 5.442 mm | 5.028 mrad | 56.105 kJ |
| 1.25 ms | 42.523 mm | 8.224 mm | 2.971 mrad | 56.158 kJ |

Phases and shared wave samples agree exactly. The surge self-change grows
rather than shrinks on the second refinement; this 30 s source prefix has
not established step convergence. Internal PTO work changes by at most
0.142% per halving. The previously independent public Python
`WEC.run(..., orifice_force_path="published_owc")` trajectory has 30 s surge
3.092634 m and PTO work 57.419 kJ. Against the finest local source these
differ by 52.719 mm and 2.244%, respectively. The Python trajectory uses
only the pinned HDF5 and phases, not MATLAB motion or force histories.

With the wave/radiation data grid fixed at 5 ms, reducing only MaxStep from
2.5 to 1.25 ms moves source surge by up to `33.996` mm over 30 s. At the
same 1.25 ms MaxStep, changing `simu.dt` from 5 to 1.25 ms leaves only
`4.829` mm surge difference. The MaxStep-only arm removes `98.666%` of the
squared surge-history gap between the uniform 2.5 and 1.25 ms arms, but only
`60.3%` of the analogous heave gap. Its internal PTO work is `56.142` kJ,
within `0.029%` of the uniform 1.25 ms arm. The solver maximum step is the
leading cause of this 30 s source surge shift; the other channels do not
share that attribution.

A separate local 5 ms source copy set only the existing flexible-state and
joint `To Workspace` blocks' `SampleTime` to `-1`, recording 76,945 strictly
increasing internal times over 30 s. Interpolated surge position, heave
speed, and flexible speed reproduce the archived 5 ms output within
`3.2e-14` in their units. Integrating the internally logged heave and
flexible hydrodynamic accelerations over those times misses their saved
speeds by at most `0.00401` and `0.01273` m/s. Coarsening those same signals
to the 5 ms output grid raises the maximum integration mismatch to `0.1513`
and `1.237` m/s. The internally logged flex acceleration is hydrodynamic
force divided by the pinned 91.61377 kg effective mass; the source's separate
reported acceleration also includes piston force. This directly identifies
output-grid aliasing in the apparent heave/flexible force imbalance. It does
not explain the 30 s Python/source motion and work gaps or establish source
step convergence. The full 130 s published case remains qualified, with no
badge increment or change to Python's physical default.
