# OWC piston-force path: corrected-source physical comparison

The pinned Applications `OWC/OrificeModel` sends hydrodynamic excitation to
its flexible `GBM` state-space input, but sends excitation **plus piston
force** to a separate reported-acceleration sum. The state-space matrix
already contains flexible stiffness and damping. On a disposable copy of the
published model, [this MATLAB script](../../../matlab/owc_orifice_correct_piston_path.m)
adds the existing piston signal to the state-space force input. It leaves the
reported acceleration route, rigid-body coupling, orifice law, HDF5, and wave
settings alone. This is a derived physical model, not the published application.

The source checkouts are WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`
and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The original
model and cleaned HDF5 SHA-256 values are
`53cdaf4830d7ac9340a17a8206b525b6c2b405ac18d8cc07912541ea2edec3cc`
and `fe702ecf99f10d4cd0bba5d858b588dd37709cd45a5c3727356849e8dc40d676`.
The locally saved corrected model is
`5df38eccf4a04aac709710b279d4bec2ee58a0bd129435e8d2910a578d05cbea`.
Simulink rewrites its archive on save, so the patch script checks the force
connection rather than requiring an identical output-file hash. The prepared
input-file SHA-256 is
`664b28b9d336b9d33d3245291f9a26471e7d64f1dd7f7498eac7f1f3cbf461f6`;
its only functional changes from the published input disable Explorer and
post-run plotting and remove unused Control System Toolbox declarations.
As a separate reproducibility check, the committed patch script was applied
to a fresh copy restored from the original model hash above. With the same
0.0003125 s MaxStep and a 15 s end time, its 3,001 body, flexible-mode, and
wave rows and 500 wave components match the earlier corrected 15 s run
element for element (maximum absolute difference zero in each CSV).

Two local R2025b Update 7 runs on Mali kept the published 130 s duration,
`simu.dt=dtOut=0.005` s, PM seed 128, 10 s ramp, 15 s radiation memory, and
`ode23t`. After initialization, only `MaxStep` differed: `0.0003125` and
`0.00015625` s. The [refined-step runner](../../../matlab/owc_orifice_refined_maxstep.m)
accepts either value as `auditMaxStep` on a patched case copy. The arms logged
416,001 and 832,001 strictly increasing internal orifice samples. Their 500
PM components and all 26,001 wave samples match exactly, including the
unmodified-model source fixtures. The first 15 s of the coarser full run
reproduces the separate corrected 15 s run byte-for-byte for body, flexible
mode, wave, and component data.

Before the full runs, source self-stability and independent Python/source
motion limits were set at 20 mm surge, 10 mm heave, 1 mrad unwrapped pitch,
and 2 mm flexible displacement. The source PTO-work change had to be at most
0.5%, and Python/source PTO-work difference at most 1%. Python independently
advances the public physical `WEC.run` route from the HDF5 and saved PM phases;
no MATLAB body state or force is an input.

| Full 130 s quantity | MATLAB step self-change | Python vs finer MATLAB | Limit |
| --- | ---: | ---: | ---: |
| Surge peak | 0.418 mm | 2.246 mm | 20 mm |
| Heave peak | 0.114 mm | 0.116 mm | 10 mm |
| Unwrapped pitch peak | 0.371 mrad | 0.359 mrad | 1 mrad |
| Flexible displacement peak | 0.020 mm | 0.022 mm | 2 mm |
| Internally integrated PTO work | 863.81 to 864.42 J; 0.071% | 865.04 J; 0.072% | 0.5% source; 1% Python |

All declared screens pass. The [compact source fixture](source.npz) retains
both full motion histories, the common wave and components, internal work,
solver settings, and Mach diagnostics (SHA-256
`eca4075bb88952931f50200ceeb0ff10ac9396afe539a47656c9be8d5818914e`).
Its coarser-arm body/orifice CSV SHA-256 values are
`d996e87f7166469a66a3e61f99c2ac2682ca83b4e0df0bd64e70d261ebdcde93`
and `c920ccb7deed2332565a0fe533c6b4c590768ebf9ed9d7e6975e0ccf354d1cfd`;
the finer-arm values are
`a73e693d7c2fb1f7ef14201b604ad16b10ddb8a5bed5f24cd2fdc53f0fa88f8e`
and `d070bec5b5dc8004a20823b2a988bd094e0f0f5628ba17bd56c1c530810f112c`.

At the **same** 0.0003125 s MaxStep, the unmodified and corrected MATLAB
models differ by up to 1.948 m heave and 0.281 m flexible displacement.
Their internal PTO work is 331,564 J versus 864 J. In the finest unmodified
run, the source's Mach warning flag is set on 62.2% of internal evaluations
and reaches Mach 1.113; in the corrected run it never rises and maximum Mach
is 0.155. That flag does not switch force laws. The one-connection edit
removes the large physical-path discrepancy; it does **not** establish parity
with the original published Simulink model. `OWC/OrificeModel` remains
`qualified_motion`, and strict badge counts do not change.
