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

The finer source output still contains rapid heave/flexible acceleration
variation; differentiating its 2.5 and 1.25 ms saved heave velocity does not
reproduce the acceleration output at those sampled times. That observation
does not identify a model defect: the force and acceleration signals can vary
between saved samples. The audit therefore leaves the full 130 s published
case qualified, with no badge increment or change to Python's physical
default. Its next useful comparison is an internal-step flexible-state/force
trace, rather than another unchanged output-grid source run.
