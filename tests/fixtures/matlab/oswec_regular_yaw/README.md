# Published OSWEC regular selected-heading yaw

`source.npz` retains all 60,001 samples of the published 600 s
`Variable_Hydro/Passive_Yaw` regular-wave case from the successful
[R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38026839689)
(run commit `a0cad6aa77d34de7b0c7e12830d3d6f39502362e`). It stores the
unrounded source time, flap yaw and speed, six excitation components, fixed
base position and velocity, PTO force, speed, and power, wave elevation, and
selected heading. Its SHA-256 is
`cd5a34dd12bdbadac653c57d7ce5f351c1f8c1b8798e09f35551c368ab5404af`.
Regenerate it from the five `OSWEC_VARIABLE_YAW_PUBLISHED` CSV files in that
run's `matlab-OSWEC_VARIABLE_YAW_PUBLISHED-baseline` artifact using the column
selection in `tests/test_oswec_regular_yaw_portable.py`.
The source body 1, body 2, PTO, wave, and selected-heading CSV SHA-256 values
are, respectively,
`a51adb66f9734bf2925aef209acdf1e03e235d831ad2bb2b593f768221cbc0fa`,
`2c0d0b2d849d6c8e9998c716bc5a0df36d4c73ac2c23ede8d94294aa658c15a8`,
`6714abb00d5061c4b406ff6c8fe5ceeb63f32ac504ac8c228b8de266e5921221`,
`35f1458959fb4a497a86c40bf2566b24f21cfd48e2ccb950f32a99ee6d90d51e`,
and `4bf2ec6d4e41fbf07da1137a688922eb4f900934e8c8467e3d2bc60a3a7eb506`.

The source HDF5 SHA-256 is
`026c9116294dcf88aff57e000a42b4b3fcfc4a3b5a79bf4269a32fc0956b1195`.
The existing bundled `oswec_multi_wave/oswec.h5` has all the same 68 dataset
paths. Every dataset except four impulse-response arrays is exactly equal;
the largest absolute difference among those four is `4.01e-11`. Using that
bundled HDF5, the independently advanced public `WEC.run` stays within the
original published-case limits over all 600 s. Both Python and MATLAB choose
the same 374 heading-change samples, and local maximum yaw, speed, and PTO
torque differences are `6.53e-15` rad, `3.68e-16` rad/s, and `1.32e-10` N m.
The relative PTO-work difference is `7.30e-14`.

Ordinary Python CI uses this NPZ and the bundled HDF5. The fresh MATLAB
workflow uses its newly generated CSVs and HDF5 for the same full-duration
gates. This strengthens an already counted published OSWEC case.
