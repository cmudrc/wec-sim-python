# Derived ellipsoid ode45 heave mass correction

This is an on-demand investigation of [issue #599](https://github.com/cmudrc/wec-sim-python/issues/599), not a replacement for either published `Nonlinear_Hydro/ode45` input. The source is WEC-Sim core `0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The [source script](../../../matlab/ellipsoid_ode45_mass_correction.m) copies each application model and input into a disposable sibling directory. It does not write to the pinned WEC-Sim or Applications repositories.

The core `bodyClass.adjustMassMatrix` shifts `D = adjMassFactor × trace(A_xyz)` into the rigid translational mass and leaves `A_zz − D` in the acceleration force sent through a 100 ns Transport Delay. The copied input sets `adjMassFactor = A_zz / trace(A_xyz)` from the same HDF5 and wave mode. For these **single active heave-coordinate** cases, the corrected model therefore has effective heave mass `m + A_zz` and zero active delayed coefficient. It retains the linked core block and the original wave, geometry, damping, PTO, and 150 s duration. The scalar shift also changes inactive surge and sway masses, so this is not a general multidimensional core fix.

The [R2025b reproduction run](https://github.com/cmudrc/wec-sim-python/actions/runs/40102870958) at Python-repository head `963f97442da8258578c1fde94c9f6d3bcb78050e` completed both corrected full-duration source arms for each wave case. The first uses the published 0.05 s solver, output, nonlinear-force, and CIC clocks. The second changes solver, nonlinear-force, and CIC clocks to 0.005 s while retaining the 0.05 s output grid. The original published input is unchanged. Both arms have 3,001 output samples and identical wave elevation at shared times. The initialized heave delayed coefficient is zero, the actual applied delayed heave force stays below `6.0e-22` N, and all inactive coordinate accelerations stay below `1.3e-21` in their respective units.

| Maximum over 150 s | Regular | RegularCIC |
| --- | ---: | ---: |
| Derived mass factor | 0.655990137167 | 0.783566018956 |
| Corrected source 0.05→0.005 s heave change | 21.268 mm | 21.414 mm |
| Unchanged source 0.05→0.005 s heave change | 17.090 mm | 17.919 mm |
| Corrected 0.05 s versus independent Python heave | 23.525 mm | 23.770 mm |
| Corrected 0.05 s versus Python speed | 29.450 mm/s | 31.858 mm/s |
| Corrected 0.05 s versus Python PTO force | 35.340 kN | 38.230 kN |
| Corrected 0.05 s versus Python absorbed power | 28.745 kW | 30.983 kW |
| Corrected 0.005 s versus independent Python heave | 2.265 mm | 2.378 mm |
| Corrected 0.005 s versus Python speed | 2.886 mm/s | 2.995 mm/s |
| Corrected 0.005 s versus Python PTO force | 3.464 kN | 3.594 kN |
| Corrected 0.005 s versus Python absorbed power | 2.804 kW | 2.983 kW |
| Corrected 0.005 s versus unchanged 0.005 s heave | 0.651 mm | 0.531 mm |

The [portable Python gate](../../../test_ellipsoid_ode45_mass_correction.py) advances motion and PTO independently from the archived source-generated `ellipsoid.h5`, the byte-identical pinned STL, and the published wave settings. Its refined 0.005 s comparisons pass the existing pointwise physical limits: 3/3.5 mm heave, 4/4.5 mm/s speed, 5/6 kN force, and 4/5 kW power for Regular/RegularCIC. Strict expected failures apply those same physical limits to the corrected published clock and retain the source-step stability limits on the corrected copy; an unexpected pass fails the test. The separate gate against each **unchanged original** published input remains unchanged and unpaired. The corrected published-clock heave misses both the 3/3.5 mm physical limits and the broader 22/23 mm numerical solver envelope used by the unchanged-input test. The corrected source's 21 mm self-change exceeds the established 3/3.75 mm heave stability screens. Removing the heave added-mass feedback alone therefore does not close the published-clock gap. Case status and badges remain unchanged.

`regular.npz` and `regular_cic.npz` archive the corrected source's time, heave, speed, PTO force/speed, wave, applied delayed heave force, inactive speeds/accelerations, and mass identity for both clocks. They were packed from the two uploaded MATLAB artifacts with [`tools/pack_ellipsoid_ode45_mass_correction.py`](../../../../tools/pack_ellipsoid_ode45_mass_correction.py); that script checks their full output grids before packing. The 882 KiB `ellipsoid.h5` here is from the Regular source artifact (SHA-256 `0092f6dd4a00c1ce6c6836df1a06bb3047ffa0e3f32038d96bcad4a408b5664c`). The two newly generated HDF5 files have different binary hashes but all 38 HDF5 datasets agree numerically, treating matching NaNs as equal. The older bundled test HDF5 has some additional datasets and minor differences in displaced volume and stored impulse-response arrays, so this portable gate uses the archived source HDF5. The STL in Applications and the bundled test geometry are byte-identical (SHA-256 `a225213457fb9273a327a141273b9d68266e650cf973a98f84e84e8607cbb7ec`).

The [workflow](../../../../.github/workflows/ellipsoid-ode45-mass-correction.yml) regenerates source evidence on request or when its diagnostic code and fixtures change, and checks fresh traces against the archived fixtures. The compact fixtures make the independent refined physical gate and both strict expected original-clock failures run in ordinary Python CI without MATLAB.
