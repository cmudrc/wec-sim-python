# MATLAB WEC-Sim comparison

This fork started from THREDgroup/WEC-Sim-Python `c04ef427` (October 2021).
The current MATLAB reference reviewed for this baseline is
[WEC-Sim/WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`](https://github.com/WEC-Sim/WEC-Sim/tree/0753b2e47f2457c078751dcfe5d251d1767b80ab).
The original Python repository includes MATLAB-generated fixtures, but there
is no MATLAB or Octave executable on the local machine. The `MATLAB reference
parity` GitHub Actions workflow checks out the pinned MATLAB revision, runs
its `waveClass` under MATLAB R2025b, and compares the generated outputs with
the production Python code. The live wave comparison passed on 6 October
2026. The workflow saves the MATLAB outputs as an artifact.

| Behavior | Reference | Baseline result |
| --- | --- | --- |
| Regular-wave elevation, ramp, and three gauge positions | Original MATLAB-generated `regular_1_test` files | Production `WaveClass` agrees within `1e-12` absolute error. |
| Finite-depth wave number | Current MATLAB `calcWaveNumber.m` dispersion relation | Python residual is checked to relative tolerance `1e-13`. |
| Finite-depth regular-wave power | Current MATLAB `waveClass.m` group-velocity equation | Corrected denominator to `sinh(2kh)`; tested at 25 m depth. |
| Current PM and JONSWAP spectra and seeded elevation | Pinned MATLAB `waveClass` with Traditional, EqualEnergy, and narrowed EqualEnergy 64-frequency discretizations, 8 s peak period, significant heights of 2.5 m and 4 m, and its Threefry phase output; [live R2025b gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37780011107) | Production `WaveClass.waveSetup` uses the current height-dependent PM and JONSWAP equations, including inferred JONSWAP `gamma`. Its explicit `phaseGenerator="matlab"` now independently generates all six saved phase matrices from `phaseSeed=1` within `1e-13` rad, including three incident PM directions. All six live wave comparisons pass locally without a phase file. In the equal-energy pairs, the maximum difference across frequency, width, and spectrum entries is `5.4e-15`; deep-water power differs by at most `2.2e-11` W/m, and origin and three-marker elevation by at most `7.4e-15` m. NumPy remains the default phase generator. |
| Narrowed PM and JONSWAP frequency ranges | Same pinned MATLAB `waveClass`, with `bem.range=[0.5,1.5]` rad/s inside a `[0.4,2.0]` rad/s BEM interval; [passing R2025b gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37780011107) | The public `PMWave` and `JONSWAPWave` settings now pass an optional range through the hydrodynamic `WEC` runner. The production component builders' selected frequencies, widths, and spectral values differ from MATLAB by at most `5.8e-15`; replayed origin elevation differs by at most `1.9e-14` m. The live gate also checks power and three marker elevations through `WaveClass`. Out-of-BEM range endpoints clamp to HDF5 limits as in MATLAB; hydrodynamic bodies retain HDF5 water depth. This validates wave preprocessing and API propagation, not a new coupled-motion case. |
| Traditional PM and JONSWAP production components | Pinned R2025b `waveClass` with 64 uniform bins, saved phases, and origin elevation from the live wave gate above | The production irregular-wave builders now reproduce both source frequency grids and spectra within `2e-12` relative and `1e-14` absolute tolerances, and origin elevation within `2e-12` m. The public `PMWave` and `JONSWAPWave` use 1000 bins by default when `discretization="traditional"`; the RM3 floating-joint runner accepts a single zero-heading irregular sea. This is wave preprocessing and a short local coupled-run check; the published RM3 MoorDyn visualization trajectory remains unpaired. |
| Irregular Bretschneider equal-energy binning | Original MATLAB-era test constants | Production `WaveClass` agrees after NumPy/SciPy API updates. This is historical compatibility only: current MATLAB WEC-Sim rejects the `BS` option. |
| Wave-surface grid for no-wave, regular, and directional irregular waves | Current MATLAB `waveClass.waveElevationGrid` equations | Implemented and checked for regular and two-direction irregular cases. |
| Published RM3 and OSWEC `Wave_Markers` elevations | Pinned MATLAB Applications inputs with 25 world XY markers, 400 s regular waves, and focused R2025b runs for [RM3](https://github.com/cmudrc/wec-sim-python/actions/runs/37934252549) and [OSWEC](https://github.com/cmudrc/wec-sim-python/actions/runs/37932186068) | Public `RegularWave.elevation_at` independently samples the source marker locations from wave height, period, direction, BEM water depth, and 100 s ramp. Against all 4,001 times and 25 markers per case, maximum elevation errors are `3.76e-11` m (RM3, 200 m depth) and `7.12e-14` m (OSWEC, 10.9 m depth), below a `1e-10` m gate. Source origin elevation and marker ordering are checked separately. Marker glyph style, size, color, and visualization files remain unpaired. |
| RM3 HDF5 hydrodynamic input | Original `rm3.h5` | Both bodies load, including their names and water depth, under NumPy 2. |
| Current RM3 and OSWEC HDF5 hydrodynamic inputs | Pinned current MATLAB core examples | Both bodies in each model load and complete regularCIC or directional irregular preprocessing with finite restoring, added-mass, excitation, and radiation arrays. This is input compatibility, not motion parity. |
| RM3 regular-wave force preprocessing | Original MATLAB-generated `body_1_test` constants and files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF; the RM3 runner now completes preprocessing. |
| RM3 irregular-wave force preprocessing | Original MATLAB-generated `body_2_test` files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF. |
| RM3 body interaction force preprocessing | Original MATLAB-generated `body_4_test` through `body_9_test` files | All six regular/regularCIC variants, including body-to-body coupling on/off and state-space radiation on/off, agree for both RM3 bodies on restoring stiffness, added mass, excitation, radiation IRF, and state-space matrices where present. |
| OSWEC directional irregular force preprocessing | Original MATLAB-generated `body_3_test` files | Production `BodyClass` agrees for restoring stiffness, added mass, radiation IRF, and three-direction excitation after replacing removed SciPy `interp2d`. |
| Sphere `noWaveCIC` heave free decay (0 m, 1 m, 1 m with Morison elements, 3 m, 5 m) | Current MATLAB WEC-Sim and MATLAB-generated Sphere HDF5 | The focused Python linear heave solver agrees over 40 s to maximum differences of 0.45 mm position, 0.64 mm/s velocity, and 363 N total force in the 5 m case. The Morison element has only x-direction coefficients in the published 1 m case, so it does not affect heave. |
| Sphere sampled-elevation heave | Pinned `Free_Decay/0m` Sphere model with only its no-wave input replaced by a reproducible two-frequency `elevationImport` record; [fresh R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37785332545) | The public `WEC.run(ImportedElevationWave(...))` path reuses the production imported-elevation convolution and runs all 4,001 samples over 40 s. With the source body block's second force ramp explicitly selected, maximum MATLAB/Python differences are `6.7e-15` m elevation, `3.7e-9` N or N m across six excitation components, 0.0243 mm heave, and 0.0323 mm/s speed. Paired gates are `1e-12` m, `1e-6` N or N m, 0.1 mm, and 0.1 mm/s. This derived case validates one-body heave dynamics with a sampled sea; the published two-body RM3 MooringMatrix imported-elevation motion remains separately paired through its floating-joint solver. |
| Fixed monopile Cartesian Morison force | Pinned MATLAB Applications `Morison_Element/morisonElement` and R2025b run, with no HDF5 bodies | The public `WEC.fixed_body` and `WEC.morison_element` configuration replays the published 400 s, three-heading PM sea and its six-component stationary-body force using `PMWave(seed=5, phase_generator="matlab")` without saved phases. All 1,500 phase values match the source within `1e-13` rad. All 500 equal-energy frequencies, spectral amplitudes, widths, and finite-depth wavenumbers agree within `1e-12`; elevation agrees within `1.9e-13` m. On all 40,001 samples, the largest force/moment difference is `7.3e-6` N m against a `1e-3` paired gate, and both body positions and velocities agree exactly. The source logs the negative of physical Morison force. Its `irregWaveMorison.m` uses the first random-phase column for every force heading even though wave elevation uses heading-specific phases; Python's explicit `phase_mode="matlab_shared"` reproduces that source behavior while the default uses heading-specific phases. This validates the published fixed Cartesian element, not moving elements or normal/tangential coefficient mode. |
| Fixed hydrodynamic monopile and tower | Pinned MATLAB Applications `Morison_Element/monopile`, generated `monopile.h5`, and [fresh R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37858193863) | The public `run_fixed_hydro_monopile` synthesizes the published 500-bin PM sea, six-component hydrodynamic excitation, two stationary-body force traces, and two fixed-joint reactions. On all 40,001 samples over 400 s, maximum elevation error is `2.1e-13` m, excitation and hydro-body total-force error `5.9e-7` N or N m, and joint-reaction error `2.6e-7` N or N m, against paired gates of `1e-10` m and `1e-5` N or N m. Both positions and velocities agree exactly. Source radiation, added mass, Morison force, damping, and acceleration are zero in this stationary case. The tower's weight and joint reactions close static force and moment balance. This verifies the published fixed hydrodynamic case, not a moving monopile or general constraint solver. |
| MBARI `Cable` coupled three-body regular-wave case | Pinned Applications `Cable`, generated `mbari.h5`, and [fresh R2025b force and joint export](https://github.com/cmudrc/wec-sim-python/actions/runs/37859525482) | The public `run_mbari_cable` independently advances both hydrodynamic bodies and the nonhydrodynamic cylinder with the published spherical joint, axial spring/damper, and drag on both cable endpoint bodies. Over all 12,001 samples and 120 s, maximum body-position error is 4.3 cm, pitch error 0.0084 rad, and body-speed error 0.422 m/s at a stiff cable snap. Cable force differs by 3.03 kN RMS and 39.9 kN at its worst instant against a 200 kN source peak; signed force-impulse error is 234 N s. The paired gates require positions below 5 cm, active speeds below channel-specific limits up to 0.45 m/s, cable force below 45 kN pointwise and 3.5 kN RMS, and impulse error below 500 N s. The source radiation, hydrostatic, and viscous body force laws reconstruct on saved states within `1e-7` N or N m; the source HDF5 added-mass and damping matrices agree with Python within `1e-8`. Given the saved body poses, rates, and accelerations, cable axial force, endpoint drag, and each 1 kg endpoint inertia reconstruct the source lower/upper vertical joint loads within 2/7 N. The independent runner retains physical implicit added mass and omits the two 1 kg endpoint inertias; the pinned Simscape source shifts added mass into rigid bodies with delayed acceleration feedback. This is planar published-case trajectory parity, not a general six-DOF cable or constraint solver. |
| WaveBot `Load_Mitigating_Controls/CalcImpedance` three-DOF system identification | Pinned Applications multisine A, WAMIT/BEMIO `waveBotBuoy.h5`, and a fresh R2025b run of the published 600 s no-wave case | `run_wavebot_impedance` independently integrates surge, heave, and pitch with 10 s radiation memory, implicit added mass, mooring stiffness/damping, quadratic drag, and the published multisine force. Source-convention motion differs from MATLAB by at most 0.13 mm surge, 0.20 mm heave, and 0.00028 rad pitch over all 60,001 samples; speed differences stay below 1.7 mm/s. Source force diagnostics reconstruct actuation, mooring, quadratic drag, and linear damping; the independent mooring-force error stays below 3.1 N. The published input uses `linearDamping(1:2:5)`, which fills the first matrix column in MATLAB and couples heave/pitch damping to surge speed. Its pitch command is applied with a negative sign. `source_linear_damping=True` reproduces that source input; the default uses a dissipative diagonal matrix and differs from the MATLAB heave/pitch trajectories by up to 17.7 mm/0.0079 rad. The source damping matrix delivers positive instantaneous power during 42% of samples, although net work is dissipative over the run. This is the open-loop system-identification case; adaptive `ControlTests` is paired below. |
| WaveBot `Load_Mitigating_Controls/ControlTests` adaptive power controller | Pinned Applications `ControlTests` input, `Z_included.mat` and `f_vec.mat`, WAMIT/BEMIO `waveBotBuoy.h5`, and [fresh 100 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37863811457) | `run_wavebot_power_control` independently synthesizes the 500-bin PM sea, integrates surge/heave/pitch with the published mooring, radiation memory, drag, and PI actuation, samples velocity and delayed command at 0.25 s, estimates complex excitation with a 1024-point Hamming FFT, and updates six gains every 8 s using the source power-only objective. Replaying MATLAB's saved signals reproduces all 13 complex excitation estimates within `2e-12` and reconstructs the PI force within `1e-8` N or N m. In the independent 10,001-sample closed loop, maximum source-convention differences are 0.093/0.228 mm surge/heave, 0.000228 rad pitch, 1.25 mm/s speed, 1.10 N or N m actuator command, 33.3 gain units, and 0.060 N or N m complex excitation estimate. The first optimizer update has a shallow, nonunique surge/pitch objective: SciPy and MATLAB can return different gains from the same estimate while attaining near-equal cost. The source input's indexed linear damping is opt-in (`source_linear_damping=True`); the Python default remains dissipative and will follow a different trajectory. The supplied impedance has bins with negative real damping, so this is numerical parity for the published controller, not a passivity claim or a general controller API. |
| WaveStar `WECCCOMP/WECCCOMP` unforced base case | Pinned Applications base input, BEMIO-generated `wavestar.h5`, and [fresh 141.2 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37866626954) | The public `run_wavestar_published` independently integrates the float and arm around their common revolute pivot, body-local B-to-C PTO stroke, 500-bin JONSWAP excitation, implicit infinite-frequency added mass, restoring and linear damping, and the published 92-state active-input radiation fit. The five-body source has a fixed frame and an unforced PTO, so those elements reduce to one moving coordinate in this case. Python reproduces frequency bins within `5e-14` rad/s, wave elevation within `9e-15` m, and six excitation components within `1.2e-12` N or N m. Across all 14,121 paired samples, maximum float position errors are 3.95 µm surge, 7.60 µm heave, and `1.57e-5` rad pitch; pitch-speed error is `1.58e-4` rad/s, PTO stroke and speed errors are 3.11 µm and 0.032 mm/s, and six-component radiation error is below 0.001 N or N m. Replaying saved MATLAB velocities through the same fit reconstructs logged radiation within 0.003 N; source force balance closes to numerical precision. The fitted active-joint damping stays positive at sampled frequencies but differs substantially in magnitude from the BEM table, so this is numerical parity for the published fit rather than a general passivity or physical-accuracy claim. The pinned input names lowercase `.stl` files while Git stores uppercase `.STL`; the MATLAB baseline corrects only those filenames in its temporary checkout. Fault and nonlinear predictive control variants remain open. |
| Fixed monopile with uniform, power-law, and linear current | Three derived 20 s, one-heading cases from pinned `Morison_Element/morisonElement`; [R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37820035413) | `PMWave(current=Current(.8, 45, profile, 30))` adds a ramped horizontal current to the fixed element's fluid velocity. Each profile has a full 2,001-sample MATLAB wave and six-component force time-series pair. The largest force/moment errors are `1.23e-6`, `1.23e-6`, and `1.16e-6` N or N m for uniform, one-seventh power, and linear profiles, respectively, against a `1e-3` gate. Wave elevation passes a `1e-10` m gate for each. Active current with multiple wave headings is rejected because the pinned MATLAB force function counts it once per heading. Regular-wave current on a moving axial element is supported separately below; its coupled trajectory is unpaired. |
| Moving Cartesian Morison source force in regular waves | Pinned `regWaveMorison.m` option 1, called at eight prescribed six-DOF states during the [MORISON_FIXED MATLAB baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37768083149) | The separate `regular_morison_source_force` diagnostic matches all six MATLAB force and moment components within `5.5e-12` N or N m, including nonzero body velocity and acceleration, angular motion, wave ramp, and an emerged zero-force state. These are prescribed-state force checks, **not coupled moving-body trajectory parity**. The pinned source rotation matrix is nonorthogonal at the tested nonzero attitudes and its angular kinematics use the unrotated local point; the diagnostic reproduces those source rules without changing the default WEC dynamics. General moving Morison layouts remain unsupported in the public `WEC` runner. |
| Moving normal/tangential Morison source force in regular waves | Pinned `regWaveMorison.m` with `bodyMorison=2` at eight prescribed moving six-DOF states; [R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38107157315) | The separate `regular_morison_source_force(..., body_morison=2, element_axes=...)` diagnostic projects fluid and body velocity and acceleration onto a nonunit body-local element axis, then applies the two drag and added-mass coefficients. Against all six source force and moment channels, the local maximum difference is `4.8e-12` N or N m (gate `1e-8`). The pinned source retains its nonorthogonal rotation convention. This is a prescribed-state source-force comparison, **not** a coupled six-DOF trajectory or support for normal/tangential elements in the public `WEC` runner. |
| Moving Morison source force with regular-wave current | The pinned `regWaveMorison.m` option-1 function was evaluated at eight prescribed six-DOF states for each of uniform, one-seventh power, and linear currents in a [fresh R2025b baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/38021058710). | The separate `regular_morison_source_force` diagnostic reproduces all six force/moment channels within `5.5e-12` N or N m for every profile, against a `1e-9` gate. The public `RegularWave(current=Current(...))` and three-coordinate moving axial Morison runner add ramped current through relative velocity; a finite-depth public run and an analytical tilted-element drag check pass locally. This is a prescribed-state source-force pair and a public-path check, **not** coupled current-trajectory parity. The public force uses a proper pitch rotation while the diagnostic retains the pinned source rotation convention. |
| Sphere moving Morison with a power-law current, coupled diagnostic | A [fresh 40 s R2025b derived run](https://github.com/cmudrc/wec-sim-python/actions/runs/38022837805) keeps the pinned `Free_Decay/1m-ME` Sphere joint, 1 m/8 s regular wave, and axial heave coefficients while adding a 0.8 m/s, zero-heading, 30 m-deep one-seventh-power current. The run exports all 4,001 six-DOF body states and direct `regWaveMorison.m` current-minus-no-current force on those same states. | The wave elevation matches Python within `3.8e-15` m. The separate Python source-convention diagnostic reproduces the direct MATLAB current-force increment within `3.7e-12` N or N m across all 4,001 states (gate `1e-6`). MATLAB option 1 adds exactly zero heave force from this horizontal current and at most 24.7 N surge force and 49.4 N m pitch moment. The physical axial Python law on those same MATLAB states adds up to 1.172 kN heave force; its current increment differs from the source by up to 58.1 N surge and 49.4 N m pitch moment. Relative to their respective no-current runs, the 40 s source trajectories change by at most 0.321 mm surge and 0.000640 mm heave, while independent Python changes by 4.965 mm surge and 1.581 mm heave. Direct Python-versus-MATLAB current-run position differences reach 16.96 mm surge, 1.903 mm heave, and `1.203e-4` rad pitch. This **does not establish coupled current parity**: the source Cartesian option-1 coefficient transform and Python proper-rotation axial loading give different force channels at nonzero pitch. The physical Python default is retained. |
| Sphere moving Morison heave free decay | The pinned `Free_Decay/1m-ME` application with its originally surge-only element changed to axial heave coefficients (`Cd=Ca=1`, area 100 m², volume 20 m³); [fresh R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37770321922) | The public `WEC.morison_element` couples body-local axial drag and acceleration-dependent added mass to a hydrodynamic heave body in still water. Over all 4,001 samples, independent Python motion differs by at most 0.325 mm heave and 0.952 mm/s speed. After the first second, physical Morison force differs by at most 43.2 N; the 40 s signed force-impulse difference is 198 N s. MATLAB logs its acceleration-feedback Morison force as zero at the first two samples and has a 41.2 kN pointwise difference from Python during startup; reconstructing MATLAB force from its *saved* acceleration only agrees within 14 N after 0.5 s. Python keeps implicit added mass rather than reproducing that source feedback transient. This derived case validates coupled single-body heave in still water. |
| Sphere moving Morison heave in regular waves | The same pinned `Free_Decay/1m-ME` Sphere case, with heave coefficients above and a 1 m, 8 s regular wave ramped over 10 s; [passing R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37773024911) | The public single-body heave solver evaluates fluid velocity and acceleration at the moving body-local element, uses relative velocity for drag, and includes the element's added mass in the acceleration solve. Across 4,001 samples over 40 s, wave elevation differs by at most `3.8e-15` m, heave by 2.26 mm, and speed by 2.52 mm/s. After the first second, physical Morison force differs from the logged source force by at most 0.935 kN over a roughly 65 kN peak. The source floating joint also permits surge and pitch: its surge reaches 1.42 m and pitch 0.0293 rad, whereas the Python configuration only moves in heave. On the source's actual six-DOF states, direct `regWaveMorison.m` replay matches the logged force within 10.1 N after startup; projection onto heave alone changes that source force by up to 0.973 kN. The Python heave force matches the pinned source function on those projected states within `5.9e-10` N, against a `1e-6` N gate. This is **reduced heave comparison**, not full three-DOF trajectory parity. The initial source acceleration-feedback force spike is excluded from the post-startup force gate. Other moving layouts and coupled current trajectories remain unpaired. |
| Sphere moving Morison in three-DOF regular waves | The derived 1 m/8 s `Free_Decay/1m-ME` case with its published surge/heave/pitch floating joint; [fresh passing R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37776059648) | The public `WEC` integrates a body-local axial element through all three active coordinates using a proper pitch rotation, relative fluid velocity, fluid acceleration, and an implicit positive-semidefinite added-mass matrix. Against the 4,001-sample, 40 s MATLAB trajectory, the independent Python run differs by at most 12.2 mm surge, 1.01 mm heave, and `9.8e-5` rad pitch; velocities differ by at most 1.32 mm/s surge and heave, and `4.6e-5` rad/s pitch. After the first second, force and moment maxima differ by 429 N surge, 149 N heave, and 176 N m pitch, with explicit paired gates on each channel. Evaluating the physical Python force law on the saved MATLAB states separates source rotation conventions from trajectory error: differences are at most 427 N surge, 91 N heave, and 180 N m pitch. The source-specific `regWaveMorison.m` diagnostic remains separate from the physical default. This validates the published three-coordinate Sphere case; general six-DOF moving elements and other coefficient layouts remain unsupported; current-loaded coupled trajectories remain unpaired. |
| Sphere `Controls/MPC` closed-loop heave and PTO | Pinned MATLAB Applications `Controls/MPC` at `d53d4d4`, pinned core at `0753b2e`, and [instrumented R2025b baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37743807105) | The public `run_sphere_mpc` solves the 400 s published case independently with 500 JONSWAP bins, a fourth-order radiation fit, 201-row AR excitation forecast, constrained 31-variable QP, delayed force-rate actuation, and physical heave radiation memory. Replayed phase frequencies, spectra, and widths agree within `1e-12`, wave elevation within `3e-12` m, and heave excitation within `1e-7` N. All prediction and QP matrices agree within `1e-10`; all 597 forecasts after a full history differ by at most `5e-7` N, and all 391 active QP first commands agree within 1 N/s using MATLAB's logged state and forecast. In an independent Python closed-loop run, maximum 40,001-sample differences are 1.011 mm physical heave, 1.501 mm/s speed, 17.2 N/s command rate, 16.8 N PTO force, and 2.89 kW source-signed instantaneous power. The MATLAB controller's internal plant is distinct from the physical body; its position is relative to equilibrium. The 2.5 m/8 s input infers JONSWAP `gamma=1`. The source's actual PTO force reaches 2.39 MN despite its 2 MN predicted constraint. This gate covers this single-body, heave-only application and published settings; general multi-body MPC remains unsupported. |
| Sphere `Mean_Drift` regularCIC application | Pinned MATLAB Applications `Mean_Drift` input and its BEMIO-generated Sphere HDF5 | The public Python `WEC` builder selects the control-surface mean-drift coefficient and runs the published 100 s surge, heave, and pitch case with 10 s radiation memory. The source excitation equals the first-order regular force plus amplitude-squared mean drift, each with one ramp factor; the surge drift component is 12.247 N after ramp. Paired maximum position differences are 0.043 mm surge, 0.551 mm heave, and `3.2e-7` rad pitch; velocity differences are below 0.081 mm/s surge, 1.738 mm/s heave, and `1.7e-7` rad/s pitch. The wave elevation matches within `3e-15` m, and reconstructed radiation forces differ by at most 0.24 N surge, 2.95 N heave, and `0.0002` N m pitch. The MATLAB linear model moves 16.41 m in surge by 100 s, so agreement is numerical parity, not evidence that the small-motion hydrodynamics remain physically accurate over that travel. |
| RM3 regular-wave heave, excitation, and PTO | Current MATLAB RM3 example and current RM3 HDF5 | A focused two-body linear heave solver with relative-motion PTO damping agrees locally over 400 s within 5.9 mm position and 6.1 mm/s velocity for both bodies. Heave excitation matches to less than `1e-5` N; PTO internal force differs by at most 8.2 kN over a 1.63 MN range. This reduced model does not cover surge, pitch, or full Simscape joint mechanics. |
| RM3 coupled surge, heave, pitch, and PTO | Current MATLAB RM3 example and current RM3 HDF5 | A four-coordinate two-body model uses the published pitched-slider geometry, body inertias, regular-wave forcing, and nonlinear rotation kinematics. Over 400 s, the two body surge positions differ by at most 36.8 mm, heaves by 4.18 mm, shared pitch by 0.000672 rad, and PTO force by 4.74 kN. It covers the canonical active DOFs but not general Simscape joint mechanics or other RM3 cases. |
| RM3 body-to-body Cases 1 and 2 | Pinned MATLAB Applications inputs and RM3 HDF5 generated with current BEMIO | Both published regular-wave cases run for 400 s with coupling off and on. Across the two, body positions differ by at most 36.8 mm surge, 4.57 mm heave, and 0.000948 rad pitch; PTO force differs by at most 5.38 kN over a 1.65 MN range. Turning on cross-body coupling substantially reduces the body-1 heave error against Case 2. These are reduced-model comparisons, not full Simscape mechanics. |
| RM3 body-to-body Cases 3 and 4 | Pinned MATLAB Applications `regularCIC` cases, with coupling off and on | Both the ordinary implicit-mass solver and the opt-in `simulink_delay` comparison are paired over 400 s. The ordinary solver differs by at most 67.2 mm surge, 3.64 mm heave, 0.00162 rad pitch, and 4.78 kN PTO force across both bodies and cases. The source-numerics setting reduces those maxima to 28.2 mm, 0.253 mm, 0.0000325 rad, and 0.292 kN; PTO stroke, speed, and power then differ by at most 0.235 mm, 0.244 mm/s, and 0.392 kW. Cases 5 and 6 use MATLAB's fitted radiation state-space model and are not validated by this convolution solver. |
| RM3 MooringMatrix with imported elevation | Pinned MATLAB Applications `Mooring/MooringMatrix` input, its 3,600 s `etaData.mat`, and BEMIO-generated RM3 HDF5 | The public `run_case` interface loads the MAT elevation, applies the configured joint surge spring, and runs the full 400 s comparison. The source wave agrees with Python interpolation and ramp to `8e-14` m. Production `BodyClass` uses the source's linear excitation-IRF interpolation and noncausal `same` convolution; an explicit opt-in reproduces the source body block's additional force ramp. All logged excitation-force components agree within `1e-3` N or N m. The logged mooring force equals `−100,000` N/m times its surge displacement within `1e-5` N. On all 40,001 samples, implicit-added-mass dynamics agree within 9.77 mm surge, 2.53 mm heave, 0.000605 rad pitch, 2.55 mm PTO stroke, 0.279 mm/s PTO speed, 334 N PTO force, and 1.17 kN mooring force. General mooring layouts remain unsupported. |
| RM3 MoorDyn coupled mooring | Pinned MATLAB Applications `Mooring/MoorDyn`, its 3,600 s elevation input, BEMIO RM3 HDF5, and native MoorDyn v2; [fresh 40,001-sample source baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37886952244) | The public Python `WEC.floating_joint(..., moordyn=MoorDyn(...), moordyn_point=spar.at(0,0,21.5))` configuration advances native MoorDyn from its own predicted connection motion and feeds the resulting six-component load back to the constrained two-body solver once per 0.01 s step. The 400 s paired run uses the full imported-elevation record and the source body block’s second force ramp; every saved source wave-force component matches to numerical precision. A local ARM64 run against the pinned MATLAB output differs by at most 0.785 mm spar surge, 0.0277 mm spar heave, 0.0147 mrad pitch, 116 N PTO force, 61.2 N surge connection force, 41.4 N heave connection force, 573 N m pitch moment, and 42.7 N fairlead tension. Explicit CI gates compare both bodies’ active positions and speeds, PTO force, six-component connection pose/velocity/force, and all three fairlead tensions over all 40,001 samples. The force is evaluated at a predicted pose once per step because legacy MoorDyn cannot roll back an implicit force trial. This validates the reduced RM3 floating joint and published mooring attachment, not arbitrary mooring layouts or the ParaView output. |
| RM3 End_Stops force law and refined-step trajectory | Pinned MATLAB Applications `End_Stops` input; R2025b runs at published 0.1 s and diagnostic 0.025/0.0125 s steps through the published 400 s duration | Effective stroke bounds are ±0.6 m, each spring is 100 MN/m, stop damping is zero, and the 0.0001 m transition is active. The source's logged stop force matches the Python smooth-stop law within `1e-5` N even at saved samples inside the transition. The 0.025-to-0.0125 s source stroke difference is at most 2.79 mm over 400 s; the published 0.1 s run differs from the refined source by up to 286 mm after contact. The opt-in Python adaptive solver with implicit added mass agrees with the 0.0125 s MATLAB trajectory over 400 s within 0.986 mm PTO stroke, 6.77 mm/s PTO speed, 98.2 kN PTO force, and 0.139% ordinary damper energy. Paired gates cover both active body positions and velocities, PTO motion and force, and integrated energy; a separate source convergence gate checks both refined MATLAB steps. |
| RM3 radiation force options: constant, convolution, and FIR | Pinned MATLAB Applications `Radiation_Force_Options`, BEMIO-generated RM3 HDF5, 12 s regular waves | All three published 500 s settings are checked against the Python floating-joint solver. The FIR path uses sampled kernel taps with a held radiation force during each RK4 step. Across the three cases, the largest differences are 29.6 mm surge, 2.11 mm heave, 0.00140 rad pitch, 3.25 mm PTO stroke, 1.55 mm/s PTO speed, 1.86 kN PTO force, and 1.57 kW PTO power. Direct FIR radiation-force differences remain under their paired gates. MATLAB FIR differs from MATLAB convolution by up to 345 mm surge, so the methods are not interchangeable for this case. The application's state-space trajectory remains unpaired because the pinned fit has negative low-frequency damping; its force law on prescribed motion is checked separately below. |
| RM3 state-space radiation force on prescribed motion | Pinned MATLAB Applications `Radiation_Force_Options` state-space setting and [focused R2025b gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37935700135) | `wecsim.radiation.replay_rm3_fitted_radiation` advances the HDF5 A/B/C channel fits from zero state using the saved six-component velocity histories, scales C by water density, and sets D to zero as the pinned MATLAB block does. Over 5,001 samples on both bodies, every radiation-force channel differs by less than `0.026%` of its source peak; the gate allows `0.03%` with a `0.0009` N floor for nearly stationary channels. This is a prescribed-motion force replay, not an independently coupled trajectory or evidence that the fitted model is passive. The default convolution solver is unchanged. |
| RM3 PTO extension float and spar free decays | Pinned MATLAB Applications `RM3_PTO_Extension` inputs and BEMIO-generated RM3 HDF5 | The no-wave floating-joint solver initializes the float at +5 m or the spar at −5 m, reproducing each published case's +5 m PTO stroke. Over 30 s, maximum active-body heave differences are 47.2 mm for float and 6.85 mm for spar; velocity differences are 86.1 mm/s and 1.87 mm/s. Both initial poses and the passive body's heave agree within numerical precision; PTO force and power are zero in Python and within `1.7e-8` in MATLAB output. |
| RM3 Multiple Condition Runs Option 1 physical sweep | Pinned MATLAB Applications RM3 input, expanded to eight scalar height × period × PTO-damping combinations | The Python `regularCIC` floating-joint solver agrees across all eight 400 s conditions within 75.1 mm surge, 0.387 mm heave, 0.000183 rad pitch, 0.413 mm PTO stroke, 0.231 mm/s PTO speed, 0.479 kN PTO force, and 0.487 kW PTO power. These scalar runs validate each physical condition; all three MCR drivers are separately checked below. |
| RM3 Multiple Condition Runs Options 1–3 orchestration | Pinned MATLAB Applications arrays, Option 2 Excel grid, Option 3 MAT file, and three actual MATLAB R2025b `wecSimMCR` runs | The Python loaders produce the same ordered eight-case table for all three published inputs. Each driver has a separate paired gate for every body's 4,001 position and velocity samples, every PTO stroke, speed, force, and power trace, the eight mean powers, and two 2×2 power matrices with period rows and height columns. Python mean absorbed powers differ from MATLAB by at most 0.155 kW or 0.059% across eight conditions in the MAT driver. The Option 2 Excel and Option 3 MAT driver outputs agree numerically exactly; Option 1 array and Option 3 MAT body positions differ by at most `1.9e-12` m. MATLAB's PTO power sign is opposite Python's positive absorbed-power convention. Phase-seed sweeps, multiple PTOs, and other postprocessing are not yet paired. |
| Public RM3 floating-joint configuration | Pinned MATLAB MCR Option 1 array condition 1 and published MooringMatrix application; fresh paired [MCR](https://github.com/cmudrc/wec-sim-python/actions/runs/37787926107) and [mooring](https://github.com/cmudrc/wec-sim-python/actions/runs/37787936962) runs | `WEC.floating_joint` routes a Python-configured two-body device through the already paired pitched-slider solver and returns named body, four-coordinate, and PTO histories. In the 400 s MCR case, maximum differences are 35.1 mm surge, 0.228 mm heave, `3.46e-6` rad pitch, 0.243 mm PTO stroke, 163 N PTO force, and 118 W absorbed power; explicit gates cover positions, speeds, stroke, force, and power at all 4,001 samples. A second public configuration reproduces the 400 s imported-elevation MooringMatrix run, including its initial spar offset and joint surge spring, under the existing 40,001-sample body/PTO/mooring gates. Both fresh MATLAB jobs pass. The builder's implicit added-mass default remains physical; source-compatible delayed feedback is explicit for the MCR comparison. This layout does not represent arbitrary Simscape joints or off-axis PTO endpoints. |
| RM3 imported-spectrum sea-state MCR | Pinned MATLAB Applications `RM3_MCROPT3_SeaState` and three actual MATLAB MCR runs | All three imported component tables agree to `4.9e-15`, wave elevations to `1.0e-13` m, and active excitation forces and moments to `3.8e-8` N or N m. The public `WEC.run(ImportedSpectrumWave(...))` path now replays all 4,001 incident-wave samples in each 400 s case within `1.0e-13` m; its general WEC dynamics are not claimed to reproduce the published four-coordinate RM3 motion. On MATLAB body velocities, the 60 s radiation convolution matches every logged force component within `5.4e-9` N or N m. The applied added-mass force reconstructs from the exported matrix and delayed acceleration within `7e-7` N or N m. The exact pitched-slider geometry and source-compatible mass feedback bring all three 400 s paired trajectories within 14.4 mm surge, 1.03 mm heave, and 0.000638 rad pitch; velocity differences are below 3.26 mm/s surge, 0.883 mm/s heave, and 0.000256 rad/s pitch. PTO stroke, speed, force, and absorbed power differ by at most 2.44 mm, 1.41 mm/s, 1.69 kN, and 1.19 kW. Mean absorbed power differs by at most 29 W. Explicit paired gates passed for the specialized MCR solver; PR #18 was merged. |
| RM3 PTO-Sim direct linear generator | Pinned MATLAB Applications `PTO-Sim/RM3/RM3_DD_PTO`, its BEMIO-generated RM3 HDF5, and a fresh 400 s `ode4` run at 0.0005 s | Both the focused Python two-heave runner and a public `WEC.pto(linear_generator=...)` configuration integrate the source block's d/q flux states and electrical angle with float and spar motion. The public configuration uses the existing body-local endpoint and coordinate maps. Source-only checks verify PTO relative motion, friction, load voltage, and mechanical/electrical power identities. Across all 40,001 saved samples, both paired gates require body heaves/velocities and PTO stroke/speed within 10 µm or µm/s, generator force and powers within 0.005 N or W, and all three phase currents/voltages within 0.0001 A/0.01 V. The [fresh source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37749845161) and local paired tests pass. This validates the published vertical two-body regular-wave layout, not general PTO-Sim electrical networks, radiation-memory coupling, or arbitrary attachment geometry. |
| RM3 PTO-Sim rectified hydraulic PTO | Pinned MATLAB Applications `PTO-Sim/RM3/RM3_cHydraulic_PTO`, PTO-Sim library, BEMIO RM3 HDF5, and [fresh 400 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37823137925) | A dedicated Python two-heave runner couples the published float and spar, 60 s radiation IRF, PM sea, compressible cylinder, rectifying valve, gas accumulators, hydraulic motor, generator, and PI load controller. Replayed sea elevation agrees within `1e-10` m; source-driven component equations agree to numerical precision. Across 40,001 independently integrated samples, maximum float/spar heave errors are 2.59/4.31 mm, PTO stroke 4.87 mm, PTO speed 3.55 mm/s, cylinder force 14.64 kN against a 646 kN source peak, cylinder pressures 0.335 MPa, accumulator pressure 1.46 kPa, shaft speed 0.045 rpm, current 0.023 A, and voltage 0.023 V. Valve-port flow differs by up to 0.0442 m³/s at sharp reversals (mean absolute error 0.000856 m³/s); integrated port-volume error stays below 0.000434 m³. The paired gates cover these channels. This validates the published vertical-slider case, not general hydraulic layouts or attachment geometry. The source PI law commands negative load resistance after 95.35 s for 76.2% of samples; numerical agreement does not establish passive electrical operation. |
| Sphere passive controller | Pinned MATLAB Applications `Controls/Passive (P)` and MATLAB-generated Sphere HDF5 | The Python heave model represents the published proportional controller as an 860,870 N s/m damper. Maximum differences over 400 s are 0.063 mm position, 0.065 mm/s velocity, 56 N controller force, and 32 W controller power after accounting for MATLAB's opposite force and power signs. |
| Sphere reactive PI controller | Pinned MATLAB Applications `Controls/Reactive (PI)` and MATLAB-generated Sphere HDF5 | The published active feedback maps to signed linear PTO stiffness and damping in the existing Python heave coordinate. Against fresh MATLAB output over 200 s, maximum differences are 8.45 mm heave, 5.49 mm/s velocity, 4.85 kN controller force, and 35.9 kW controller power. MATLAB's logged force and power satisfy the published gain equations to numerical precision. The 22.74 m peak-to-peak heave is the published unconstrained linear controller's response; trajectory agreement does not establish a feasible PTO stroke or physical design. |
| Sphere reactive controller with simple direct-drive PTO | Pinned MATLAB Applications `Controls/ReactiveWithPTO`, current Sphere HDF5, and a fresh 200 s ode45 run | The public Python PTO configuration couples the published PI force request to a winding L/R torque state, gear ratio, drivetrain inertia and friction, and generator current/voltage/loss outputs. On 20,001 samples, MATLAB's logged controller, shaft torque decomposition, and electrical identities agree with the source block equations to numerical precision. Its Derivative block logs shaft inertia torque from the preceding output interval's speed difference, matching that finite difference within `3e-10` N m; Python uses instantaneous acceleration in the coupled dynamics. The independent Python trajectory differs by at most 0.076 mm heave and 0.051 mm/s velocity; shaft torque by 0.350 N m and source-signed electrical power by 22.5 W. This validates one pure-heave, zero-heading regular-wave direct-drive connection. Voltage/current limits and other PTO-Sim components remain unsupported. |
| Sphere declutching controller | Pinned MATLAB Applications `Controls/Declutching` and MATLAB-generated Sphere HDF5 | A configurable body-local PTO endpoint uses the published 232,020 N s/m engaged damping and switches off for 0.8 s after velocity reversal. In 40,001 paired samples over 400 s, heave differs by at most 1.99 mm and velocity by 4.89 mm/s. Both runs have 83 disengagement events, each 0.8 s long; event times differ by at most one 0.01 s sample, affecting four force samples. Where both runs are engaged, force differs by at most 1.14 kN and power by 1.06 kW; both report zero force while disengaged. Absorbed energy differs by 10.5 kJ over a 22.4 MJ MATLAB total. A full-sample force maximum is 106 kN because one-step mode differences occur at two transitions, so force parity is gated by event timing and within-mode values. |
| Sphere latching controller | Pinned MATLAB Applications `Controls/Latching` and MATLAB-generated Sphere HDF5 | The public `LatchingControl` applies the published 49,181 N s/m normal damping and 37,308,296 N s/m latch damping for 2.4 s after velocity reversal. Both force laws dissipate energy. In 40,001 paired samples over 400 s, MATLAB and Python each have 83 latch starts; start times differ by at most 0.05 s, and each completed latch lasts 2.41 s in the sampled source implementation. Maximum heave and velocity differences are 0.379 m and 0.573 m/s; the heave extrema differ by at most 0.033 m. Absorbed energy differs by 2.43 MJ against MATLAB's 299.98 MJ (0.81%). This switched high-gain case has 553 samples in different modes, so pointwise controller force is not a useful single parity metric; the paired gate checks the exact force and power laws, event timing, motion, envelope, and integrated energy. The published linear model does not establish a feasible physical stroke or rigid latch. |
| Configured Sphere spring, damper, and PTO attachment | Derived from the pinned Sphere passive case with 50,000 N/m PTO stiffness, 100,000 N s/m native PTO damping, and the PTO moved to x = 1 m | A fresh MATLAB run and the Python `WEC` builder agree within 0.051 mm position and PTO stroke, 0.054 mm/s velocity and PTO speed, 52 N combined force, and 27 W combined power. The x offset has no effect on heave-only motion, so this validates the force settings but not attachment-location dynamics. |
| OSWEC PM equal-energy waves, directional excitation, pitch, and PTO | Current MATLAB OSWEC example and current OSWEC HDF5 | Python recreates all 500 PM equal-energy bins from the HDF5 range, then synthesizes wave elevation and six-component excitation using the saved MATLAB phase matrix. Local component differences are below `6e-15`, elevation below `2e-13` m, and excitation below `1e-7` N. The public two-body `WEC.fixed_hinge` call passed a [fresh pinned R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37795546236) with source phases: maximum flap differences over 4,001 samples were 10.0 mm surge, 1.66 mm heave, 0.00201 rad pitch, 0.00297 rad/s pitch speed, and 35.6 N m PTO torque; the fixed hydrodynamic base agreed exactly. A Python seed produces a separate reproducible realization. This model covers pitch about a fixed hinge, not a general six-DOF device. |
| OSWEC PTO-Sim hydraulic adjustable rod and fixed crank | Pinned MATLAB Applications `PTO-Sim/OSWEC/OSWEC_Hydraulic_PTO` and `OSWEC_Hydraulic_Crank_PTO`, BEMIO OSWEC HDF5, and [fresh 400 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37831809681) | Python couples the published hinged flap and 30 s radiation IRF to either rotary-to-linear linkage, compressible cylinder, rectifying valve, gas accumulators, motor, generator, and PI load. The pinned crank laws map angle to stroke and cylinder force to torque with opposite signs; source-driven cylinder pressure steps, valve flows, accumulator pressures, motor and generator steps agree to numerical precision. Across both independently integrated 40,001-sample trajectories, maximum pitch/speed errors are 0.00136 rad and 0.000722 rad/s, flap-center surge/heave errors 6.55/2.75 mm, crank torque 57.9 kN m against a 2.10 MN m source peak, cylinder force 20.0 kN, chamber pressure 0.390 MPa, accumulator pressure 5.24 kPa, shaft speed 0.0595 rpm, current 0.167 A, and voltage 0.175 V. Valve-port flow differs by up to 0.0540 m³/s at reversals; accumulated port-volume error stays below 0.00129 m³. Explicit paired gates cover body motion, torque, hydraulic states, valve flow, and electrical outputs for each layout. The adjustable-rod source commands negative load resistance after 46 s for 77.0% of its nonzero-current samples; this is a source-case limitation. General hydraulic networks and arbitrary attachment geometry remain unpaired. |
| OSWEC `Multiple_Wave_Spectra` with two independent PM seas | Pinned MATLAB Applications `Multiple_Wave_Spectra` input and BEMIO-generated OSWEC HDF5 | The public `run_case` interface sums separate 2 m, 0° and 1 m, 90° seas with saved independent phase realizations over 100 s. Combined elevation agrees within `1e-11` m and both bodies' six excitation-force components within `1e-6` N or N m. The fixed hydrodynamic base remains stationary. With the ordinary implicit added-mass default, maximum flap differences are 106 mm surge, 10.8 mm heave, and 0.0213 rad pitch; PTO torque differs by at most 7.25 kN m. The explicit source-delay comparison reduces these to 10.3 mm, 1.20 mm, 0.00207 rad, and 818 N m. Both dynamics paths have separate paired gates. The source enables passive-yaw preprocessing, but its hinge has zero yaw motion; the source's spline frequency interpolation is selected explicitly for this case. Moving-yaw dynamics and fixed-base reactions remain unsupported. |
| OSWEC `Full_Directional_Waves` imported spectrum | Pinned MATLAB Applications `Full_Directional_Waves` input and BEMIO-generated OSWEC HDF5 | The imported 47-frequency, 180-heading spectrum now generates all 8,460 source phases from `seed=1, phase_generator="matlab"` within `1e-13` rad and reproduces all 8,001 wave-elevation samples within `9e-14` m. The source force block omits heading-bin width although its wave-elevation block includes it: logged excitation is about `1/sqrt(2° in radians) = 5.352` times the physically integrated force. Python keeps heading integration as its default. Explicit `matlab_omitted` quadrature plus `spline_frequency` interpolation matches all six logged excitation components for both bodies within `4e-7` N or N m. With that source forcing and ordinary implicit added mass, the 400 s flap differs by at most 12.8 mm surge, 3.52 mm heave, 0.00258 rad pitch, and 56.6 N m PTO torque. The separate opt-in Simulink-delay comparison through public `WEC.run(FullDirectionalSpectrumWave(...))` reduces these to 1.53 mm, 0.293 mm, 0.000306 rad, and 5.54 N m; the paired gate checks the full flap, base, and PTO trajectories. Source-compatible agreement does not establish that the omitted-width force is physically correct. The fixed base remains stationary; its ground reactions are not computed. |
| OSWEC passive yaw off and on | Pinned MATLAB Applications `PassiveYawOFF` and `PassiveYawON` inputs, BEMIO-generated OSWEC HDF5, and fresh R2025b trajectories | The public `WEC` interface defines a yawing flap, fixed hydrodynamic base, and 120,000 N m s/rad rotational PTO. Both published 600 s cases use 2.5 m, 8 s waves at 10°. With passive yaw off, Python and MATLAB flap yaw agree within `4e-14` rad. With passive yaw on, excitation follows the body-relative wave heading and turns the flap toward the 10° incident heading; yaw and angular-speed differences stay below 0.00149 rad and 0.000193 rad/s. PTO torque and absorbed power differ by at most 23.2 N m and 0.128 W, accounting for MATLAB's negative absorbed-power sign. On the saved MATLAB yaw trajectory, all six reconstructed excitation components differ by at most 100 N or N m except yaw moment, which differs by 467 N m. The source holds heading coefficients until yaw changes by 0.01°, while Python interpolates continuously. The fixed base is stationary; full three-dimensional joint mechanics remain unsupported. |
| OSWEC variable-hydrodynamics passive yaw, derived 2° bank | Pinned MATLAB Applications `Variable_Hydro/Passive_Yaw` regular-wave, 600 s, 0.01 s case with its direction-bank generator resampled from the published 0.05° grid to `-40:2:40`; [fresh passing R2025b and Python run](https://github.com/cmudrc/wec-sim-python/actions/runs/38026261576) | `WEC.body(..., passive_yaw=True, yaw_heading_bank=range(-40, 41, 2))` follows the source's nearest-heading selection, single-direction `interp1(...,'spline')` frequency coefficients, and unrotated variable-hydro force convention. On all 60,001 saved MATLAB states, heading selection is identical and all six excitation components differ by at most `3.66e-7` N or N m, below the `1e-6` gate. An independent Python run from the same wave and HDF5 agrees over the full 600 s within `4.9e-15` rad flap yaw, `4.3e-16` rad/s yaw speed, and `1.9e-10` N m PTO torque on local ARM64; PTO work differs by `7.2e-14` relative. Paired gates are `1e-6` rad yaw, `1e-6` rad/s speed, `0.01` N m torque, and `1e-6` relative PTO work. The original 0.25° selected bank is paired separately below; other variable-hydrodynamics layouts remain unpaired. |
| OSWEC `Variable_Hydro/Passive_Yaw`, published 0.25° selected bank | Pinned published 600 s regular-wave `wecSimInputFile.m` and [fresh R2025b and Linux paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/38026839689) | The input's `bemDirections=-30:0.25:30` is unchanged. The pinned BEMIO interpolation generates only those selected HDF5 headings instead of all unused 0.05° files; each selected coefficient uses the same two original 10° neighbors and interpolation weight as the published generator. Python configures the same 241-heading bank and independently advances yaw, radiation, and PTO over all 60,001 samples. Heading selection matches every saved source state. All six source-path excitation components differ by at most `3.66e-7` N or N m (gate `1e-6`). On local ARM64, maximum flap-yaw, yaw-speed, and PTO-torque differences are `6.53e-15` rad, `3.68e-16` rad/s, and `1.33e-10` N m; 600 s PTO work differs by `7.45e-14` relative. Linux passed the explicit `1e-6` rad, `1e-6` rad/s, `0.01` N m, and `1e-6` relative-work gates. This pairs the published fixed-base, pure-yaw regular-wave case; general variable-hydrodynamics layouts remain unsupported. |
| OSWEC `Variable_Hydro/Passive_Yaw`, irregular PM wave, 120 s | The pinned application input with `waveFlag='irregular'`, phase seed 1, its 241 selected headings, and only `simu.endTime` changed from 600 to 120 s in a temporary copy; [fresh R2025b and Linux paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/38033071828) | Python's independently generated 500 PM frequency, amplitude, width, and phase entries differ from MATLAB by at most `5.33e-15`; wave elevation differs by `4.79e-14` m. On all 12,001 saved source yaw states, the selected heading is identical and every excitation channel differs by at most `8.06e-8` N or N m (gate `1e-6`). An independent Python 120 s run differs by at most `0.00155` rad yaw, `0.000239` rad/s speed, and `28.6` N m PTO torque; absorbed PTO work differs by `0.360%`. Paired gates are `0.002` rad, `0.0003` rad/s, `35` N m, and `0.5%` work. Both trajectories switch headings 181 times, but their selected headings first differ at 54.73 s and on 469 later samples, showing event-time sensitivity rather than a source-force mismatch. The full published run is audited separately below. |
| OSWEC `Variable_Hydro/Passive_Yaw`, published irregular PM wave, 600 s | The pinned application's unchanged 600 s wave, PTO, heading selection, and solver, with only unused 0.05° BEM files skipped by generating the same selected 0.25° coefficients; [fresh R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38035515976) | Across all 60,001 samples, Python's independently generated PM components differ by at most `5.33e-15`, wave elevation by `2.01e-13` m, and the selected heading on every saved MATLAB yaw state agrees exactly. Six excitation channels evaluated on those MATLAB states differ by at most `2.28e-7` N or N m, passing a `1e-6` source-force gate. The independent native Python run begins with close yaw agreement but differs by up to `0.1103` rad yaw, `0.01225` rad/s yaw speed, and `1.469` kN m PTO torque over 600 s. Its PTO work is `2371.33` J versus MATLAB's `2607.39` J, a `9.05%` difference. The first heading mismatch occurs at 54.73 s; MATLAB switches headings 624 times and Python 518 times. In a separate **prescribed-heading diagnostic**, Python uses only MATLAB's saved heading choices while advancing its own dynamics and evaluating its own excitation: maximum yaw, speed, and torque differences fall to `0.00868` rad, `4.06e-5` rad/s, and `4.88` N m; PTO work differs by `0.0081%`. The diagnostic gates are `0.01` rad, `1e-4` rad/s, `10` N m, and `0.02%` work. This isolates heading-event feedback as the dominant source of the long-run difference, but **native 600 s trajectory parity is not established**. The force law and default dynamics are unchanged. |
| OSWEC irregular passive yaw | Pinned MATLAB Applications `PassiveYawRegression` input, BEMIO-generated OSWEC HDF5, and fresh R2025b trajectory | The public `PMWave` replays the published 500-frequency phase realization over 250 s with 40 s radiation memory. Frequency bins, spectral amplitudes, and widths agree within `6e-15`; elevation agrees within `1e-13` m. The public `passive_yaw_threshold=1` selects the published coefficient hold and nearby BEM-heading snap. On MATLAB's saved yaw path, that opt-in force law matches all six logged excitation components within `1e-7` N or N m. The source's yaw radiation convolution reconstructs from its saved speed and pinned HDF5 within `5.4e-9` N m; added mass, hydrodynamic force balance, and rigid inertia balance close within `1.4e-8` N m. In a native Python run a heading update occurs one 0.01 s sample early at 53.45 s: MATLAB's relative heading is 10.99980° and Python's is 11.00004°, on opposite sides of the 11° update boundary. The yaw difference there is `4.3e-6` rad, and replaying MATLAB's logged force produces nearly the same pre-event drift (`4.29e-6` rad). Later event differences accumulate to 0.203 rad maximum yaw and 0.0277 rad/s speed differences over 250 s. With continuous heading interpolation, Python differs from the published trajectory by up to 0.452 rad yaw and 0.101 rad/s yaw speed. Replaying MATLAB's logged six-component force through the Python dynamics reduces these maximum differences to 0.00180 rad and 0.0000935 rad/s. The source force laws and balances are paired, but published-case trajectory parity with native excitation remains **unestablished** because the threshold amplifies small integration differences. |
| OSWEC irregular passive yaw, source-heading schedule replay | The same pinned published `PassiveYawRegression` 1° hold and three R2025b phase seeds in the [fresh three-seed paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/38024726665) | Python receives only the coefficient-heading update sequence reconstructed from each saved MATLAB yaw path; it evaluates the incident force at its own yaw and independently advances radiation, rigid motion, and the PTO over all 25,001 samples per seed. The existing source-yaw replay establishes that this heading sequence reproduces all six logged MATLAB excitation channels within `1e-6` N or N m. Across the three independent heading-schedule runs, maximum yaw and yaw-speed differences are `0.00180` rad and `9.34e-5` rad/s, and PTO damper work differs by at most `0.0313%`; explicit gates are `0.003` rad, `0.0002` rad/s, and `0.1%` work. This isolates heading-event timing as the dominant source of the native held-heading divergence. It is a **prescribed-event diagnostic**, not native full-trajectory parity. |
| OSWEC irregular passive yaw, continuous-heading control | Same pinned MATLAB Applications `PassiveYawRegression` input and three phase seeds, with only `body(1:2).yaw.threshold` changed from 1° to 0° in temporary application copies; fresh R2025b source runs for [seed 1](https://github.com/cmudrc/wec-sim-python/actions/runs/37720514831), [seed 2](https://github.com/cmudrc/wec-sim-python/actions/runs/37925418638), and [seed 3](https://github.com/cmudrc/wec-sim-python/actions/runs/37921861858) | Against the saved phases, Python advances its own 250 s yaw, radiation, and PTO dynamics. Across all three seeds, wave elevation differs by less than `1e-13` m; the six excitation components evaluated on each MATLAB yaw path differ by at most `1.15e-7` N or N m. Over all 25,001 time samples per seed, maximum native flap/PTO angle, yaw/PTO speed, torque, and source-signed power differences are `3.17e-5` rad, `7.30e-6` rad/s, `0.875` N m, and `0.0654` W. The latter is `0.0215%` of seed 2's `304.7` W peak. PTO damper work differs by at most `0.0115%`, and the fixed base remains stationary. The seed-2 source run passed MATLAB generation but its first Python check failed a seed-1 absolute `0.05` W power limit; a scale-aware `0.05%`-of-source-peak power gate plus `0.02%` net-work gate pass locally and in CI against all three saved artifacts. This control supports continuous-heading dynamics parity across these seas; it does not establish full parity for the published 1° hold. |
| OSWEC fixed nonhydrodynamic base | Pinned MATLAB Applications `Nonhydro_Body` case and its BEMIO-generated OSWEC HDF5 | The public `WEC.fixed_hinge(flap, base=wec.fixed_body(...))` call passed a [fresh pinned R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37795143695) and reports the stationary base and nonlinear flap motion about the PTO hinge. Across 4,001 samples over 400 s, maximum flap position differences are 20.8 mm surge, 7.9 mm heave, and 0.00444 rad pitch; velocity differences are 15.6 mm/s surge, 6.8 mm/s heave, and 0.00328 rad/s pitch. All six excitation-force components differ by less than 1 N, the base position and velocity agree exactly, and zero PTO torque differs only by MATLAB numerical noise below `3.3e-7` N m. The base's ground-constraint reaction forces are not calculated. |
| Ellipsoid instantaneous nonlinear hydro, `ode4/Regular` | Pinned MATLAB Applications `Nonlinear_Hydro` and BEMIO-generated ellipsoid HDF5/STL | The Python builder combines mesh buoyancy, instantaneous Froude–Krylov correction, quadratic heave drag, BEM linear excitation/radiation, and a configured PTO. On all 3,001 saved MATLAB states, mesh equilibrium mass agrees within `2e-9` kg, buoyancy within `6e-9` N, drag within `7e-11` N, and total heave excitation within `1e-6` N. The independent 150 s Python trajectory differs by at most 5.16 mm heave, 5.86 mm/s velocity, 5.16 mm PTO stroke, 7.03 kN PTO force, and 7.21 kW absorbed power. MATLAB splits the BEM added mass between Simscape mass and an applied force; Python uses the combined effective mass. This case validates the single-body, pure-heave, zero-direction constant-radiation mode. |
| Ellipsoid instantaneous nonlinear hydro, `ode4/RegularCIC` | Pinned MATLAB Applications `Nonlinear_Hydro` RegularCIC input and BEMIO-generated ellipsoid HDF5/STL | The same Python mesh model uses a 60 s radiation convolution. On all 3,001 saved MATLAB states, mesh buoyancy and quadratic drag agree within `6e-9` and `7e-11` N, linear plus nonlinear heave excitation within `4e-7` N, and radiation convolution within `3e-10` N. The independent 150 s Python trajectory differs by at most 6.39 mm heave, 7.41 mm/s velocity, 6.39 mm PTO stroke, 8.89 kN PTO force, and 8.98 kW absorbed power. The memory-step fixed-point iteration was allowed more iterations to converge near the moving waterline; no force coefficient was tuned. |
| Ellipsoid `ode45/Regular` and `ode45/RegularCIC` | Pinned MATLAB Applications `Nonlinear_Hydro` ode45 inputs and BEMIO-generated ellipsoid HDF5/STL; [full 150 s refined-step paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/38110743157) | The published 0.05 s source applies mesh buoyancy from the preceding output state. Its restoring log matches the prior-state mesh law within `6e-9` N, while the current-state mismatch exceeds 54 kN; Python retains instantaneous force evaluation. A derived comparison changes only MATLAB `simu.dt` to 0.01 and 0.005 s and independently runs Python at both steps, checking all 30,001 samples at 0.005 s. MATLAB 0.01-to-0.005 s heave changes are 1.863 mm (Regular) and 1.945 mm (RegularCIC); Python self changes are 0.0101 and 0.0107 mm. Against the finest MATLAB records, Python differs by at most 1.855 and 1.933 mm heave and 3.256 and 3.504 kN PTO force, respectively; position, speed, force, and power gates all pass. The MATLAB restoring log still reconstructs from the preceding saved state within `6.2e-9` N at both refined steps. The original 0.05 s MATLAB trajectories differ from their own 0.005 s trajectories by about 17–18 mm heave, so **published-step ode45 trajectory parity remains unestablished**. Refining `simu.dt` changes ode45 MaxStep, output, nonlinear-force, and convolution intervals together; the test does not isolate their individual effects or establish source step convergence. No Python physical default or force coefficient was adjusted. |
| Sphere variable draft and mass in heave | Pinned MATLAB Applications `Variable_Hydro/Variable_Mass`, nine BEMIO-generated draft HDF5 files, and its full 900 s ode4 run | The Python `WEC.variable_body` selects the published state every 100 s and changes equilibrium draft, rigid mass, added mass, radiation damping, restoring stiffness, and excitation together. On MATLAB's 9,001 saved states, state-selected excitation differs by at most `3e-8` N, radiation damping by `3e-10` N, restoring plus weight by `8e-9` N, and assembled hydrodynamic force by `9e-9` N. The independent Python trajectory differs by at most 2.23 mm heave, 4.14 mm/s velocity, and 827 N PTO force over all nine states. The source reaches −6.63 m by 900 s; agreement with its linear coefficients does not establish physical accuracy at that excursion. Other variable-hydrodynamic motions remain unsupported. |
| Barge generalized body modes | Pinned MATLAB Applications `Generalized_Body_Modes`, its BEMIO-generated ten-DOF barge HDF5, and a fresh 400 s R2025b run including the model's separate `Flex_out` signal | `BodyClass` pairs flexible effective mass, stiffness, damping, hydrostatic, added-mass, radiation-damping, and excitation matrices; the largest absolute matrix difference is `1.91e-6` in applied added mass. The source flexible acceleration equation and six signed force components reconstruct from those matrices, with restoring and radiation differences below `1e-6` N, and sampled added-mass feedback below `0.05` N. The public `WEC.floating_gbm` interface runs a coupled three-rigid/four-flexible-coordinate implicit-mass model in the published 2 m, 8 s regular wave. Over all 8,001 samples, maximum rigid position and velocity differences are 1.06 mm and 0.838 mm/s; flexible displacement, velocity, and acceleration differ by at most 0.043 mm, 0.034 mm/s, and `1.26e-4` m/s². This validates one single-body floating 3-DOF joint with four modes, zero-heading regular waves, and no PTO or mooring; other GBM layouts and wave/radiation options remain unsupported. |
| OWC `OrificeModel` wave, force, and early coupled motion | Pinned MATLAB Applications `OWC/OrificeModel`, BEMIO-cleaned seven-DOF HDF5, and a fresh 130 s R2025b `ode23t` [source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37799665808) | `wecsim.OrificePTO` reproduces the embedded Simulink orifice block on its 418,828 logged internal solver evaluations. Maximum force, pressure, and source-power differences are `4.8e-10` N, `2.9e-12` kPa, and `4.2e-13` kW; the Mach flag is identical. Interpolating source flow to the independent 26,001-sample flexible-mode output gives a piston-speed difference below `3.6e-13` m/s. The 500 PM bins agree within `5.2e-15`, wave elevation within `2.8e-14` m, and all six rigid plus one flexible excitation channels within `1.5e-10` N. Flexible mass/stiffness/damping, added mass, and hydrostatic matrices agree within `2.3e-13` in their units. With the [rigid-force diagnostic run](https://github.com/cmudrc/wec-sim-python/actions/runs/37804483371), radiation reconstructed from all seven MATLAB velocities differs by at most `0.0074` N in the flexible channel and `0.22` N across rigid channels (peak rigid magnitude `1,182` N). The rigid restoring and linear-damping forces reproduce the source logs within `8.6e-11` N and `1.1e-12` N, and the rigid hydrodynamic force sum closes within `1.7e-10` N. MATLAB's flexible hydrodynamic force sum closes within `1e-12` N; its effective-mass acceleration including the orifice force closes within `3.5e-8` N. The source rigid surge, heave, and pitch balances close within `1.2e-10` in their force units when using its adjusted translational mass and logged piston reaction. `WEC.floating_gbm(..., orifice=...)` now integrates the seven-DOF PM case with memory radiation, quadratic drag, and the coupled piston force. Against the saved source phases, maximum surge, heave, unwrapped-pitch, and flexible-position differences through 6 s are `0.0113` m, `0.0244` m, `0.0098` rad, and `0.0070` m. At 130 s the physical coupled path has heave and flexible-displacement differences of `1.94` m and `0.280` m. The separate opt-in `orifice_force_path="published_owc"` follows the pinned Simulink routing, which leaves the piston out of the flexible state-space input but retains it in reported acceleration and rigid added-mass feedback. Over all 26,001 samples, its maximum surge, heave, unwrapped-pitch, and flexible-position differences are `0.174` m, `0.0367` m, `0.0108` rad, and `0.0054` m; PTO reaction and absorbed-power maxima differ by `1.98` kN and `1.75` kW, and integrated absorbed work by `4.76%`. These are explicit 130 s source-convention gates, **not physical coupled-path parity**. Python reports SI watts while the source's numeric power is kW. The source first exceeds its Mach threshold at 10.13 s; it is above 0.3 for 61.4% of simulated time and peaks at Mach 1.11. Matching its incompressible law after that point does not establish physical accuracy. |
| Case-driven dynamics runner | Current MATLAB RM3 and OSWEC examples plus Sphere and RM3 Applications cases | One generalized-coordinate engine assembles rigid inertia, hydrodynamic added mass and radiation, hydrostatic restoring, excitation, and linear PTO forces. The heave, fixed-hinge, and floating-joint layouts run the paired cases above, including imported elevation, a joint surge spring for RM3 MooringMatrix, and native MoorDyn coupling for the published RM3 mooring. A fourth `linear_subspace` layout maps independent coordinates into arbitrary bodies; its RM3 two-body heave, Sphere free decay, configured Sphere PTO, and two instantaneous nonlinear-hydro cases have paired MATLAB checks. Arbitrary Simscape layouts, general moorings, and other nonlinear-hydro/application cases remain unsupported. |

The irregular variable-yaw bank is also sensitive to time-step refinement in
Python. With the published sea and bank but `dt=0.005` s, the 120 s run stops at
86.08 s because the implicit radiation-memory step alternates between the
adjacent −8.75° and −8.50° excitation headings across their −8.625° boundary.
Thirty successive trial accelerations repeat this two-cycle rather than
converging. The pinned application's `calcIndex.m` selects the nearest heading
on every evaluation; its previous-index hold rule is commented out. The
passing 120 s pair at the published `dt=0.01` s therefore does not establish a
step-converged physical trajectory. Resolving this discrete event requires a
defined within-step heading rule; no force coefficient was changed to suppress
the failure.

The published one-degree irregular passive-yaw hold has a [three-seed R2025b
audit](https://github.com/cmudrc/wec-sim-python/actions/runs/37915267886).
Each 250 s MATLAB run changes only `waves.phaseSeed` from 1 to 2 or 3;
Python uses its saved phases but advances yaw, radiation, and PTO independently.
Wave elevation differs by less than `1e-13` m for every seed, and the first
50 s yaw error remains below `4e-7` rad. Replaying the source yaw through the
Python held-heading force law reproduces all six MATLAB excitation components
within `1.13e-7` N or N m. Given the logged source force, independent Python
dynamics stay within `0.0018` rad yaw and `0.000094` rad/s speed, with PTO
damper work within `0.04%` across all three runs. With native heading feedback,
however, the maximum yaw errors are `0.203`, `0.099`, and `0.241` rad, while
250 s PTO-work differences are `0.9%`, `40.4%`, and `83.6%` for seeds 1–3.
The seed-1 work agreement therefore does **not** generalize. This is direct
evidence that the discrete heading events can amplify small trajectory
offsets; it does not establish which numerical path is physically preferable
or full published-case parity.

The seed-2 Python held-heading result is also sensitive to time step: with
the same sea phases and configuration, PTO damper work is `2,795`, `3,188`,
and `2,681` J at `dt = 0.01`, `0.005`, and `0.0025` s. For comparison, the
continuous-heading Python variant gives `4,420.58` and `4,419.92` J at the
first two steps, a `0.015%` change. Refining `dt` also refines incident-force
sampling and radiation convolution, so this experiment alone cannot assign
the entire sensitivity to the heading switch. It does show that the held
result at the published step should not be treated as a step-converged
physical estimate.

A fresh [pinned R2025b step audit](https://github.com/cmudrc/wec-sim-python/actions/runs/38094457574)
changes only `simu.dt` from the published `0.01` s to `0.005` s for phase
seeds 1 and 3. The MATLAB wave elevation is identical at common samples,
and independent Python waves differ from it by less than `1e-13` m. For
both seeds, the new published-step MATLAB yaw agrees with the earlier saved
baseline to within `5e-16` rad, confirming the repeated source run. For
seed 1, MATLAB PTO work changes from `9,068.29` to `4,735.91` J (`-47.8%`);
Python changes from `8,986.61` to `8,324.26` J (`-7.4%`). For seed 3, MATLAB
changes from `4,374.14` to `3,532.91` J (`-19.2%`); Python changes from
`8,032.48` to `3,697.78` J (`-54.0%`). The MATLAB coarse and fine yaw traces
first differ by `0.01` rad at 138.49 s and 90.48 s respectively, and differ
by as much as `0.2202` and `0.0937` rad over 250 s. Thus the published
MATLAB step itself is not a stable full-trajectory target for these held-
heading seeds. The seed-1 Python/MATLAB PTO-work match at the published step
disappears on refinement; the seed-3 work gap shrinks to `4.7%`, but its
maximum yaw difference remains `0.0953` rad. The earlier source-force and
prescribed-heading comparisons still hold. Neither implementation has a
demonstrated converged native held-heading trajectory; no Python dynamics
coefficient was adjusted to fit the coarse MATLAB output.

A [third pinned step](https://github.com/cmudrc/wec-sim-python/actions/runs/38095725988)
at `0.0025` s keeps the incident wave **exactly** equal to the `0.005` s
wave at common samples. Seed-1 MATLAB PTO work rises from `4,735.91` to
`6,246.54` J (`+31.9%`), and its maximum common-grid yaw change is `0.1216`
rad. Seed-3 work rises from `3,532.91` to `7,136.69` J (`+102.0%`), with
`0.2931` rad maximum yaw change. Independent Python held-heading work changes
from `8,324.26` to `8,231.67` J (`-1.1%`) for seed 1 and from `3,697.78`
to `4,651.48` J (`+25.8%`) for seed 3; seed-3 maximum Python self-change in
yaw is `0.2703` rad. Neither held-heading implementation has demonstrated
250 s step convergence at the finest tested step. This also means the
seed-3 `4.7%` cross-implementation work agreement at `0.005` s was
incidental, not a parity gate. The pinned
[`irregnLYaw.m`](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/functions/simulink/model/irregnLYaw.m)
describes the threshold as a way to skip coefficient interpolation and uses
zero to interpolate every step. Python's continuous-heading option is the
stable control here: across seeds 1 and 3, halving its step from `0.01` to
`0.005` s changes PTO work by at most `0.0033%` and yaw by at most
`4.5e-5` rad, consistent with its separately paired 250 s MATLAB cases.

The pinned OWC Simulink model has a force-path difference:
its flexible `GBM` state-space block integrates inverse effective mass times
the hydrodynamic force, while the piston force joins a separate reported
acceleration path. Over the first 6 s, a finite difference of the saved
flexible velocity agrees with hydrodynamic force divided by effective mass to
`0.00196` m/s² RMS; it differs from the reported acceleration that includes
the piston by `0.257` m/s² RMS. At 2 s these three values are `-0.00317`,
`-0.00299`, and `-0.15217` m/s², respectively. The rigid-heave added-mass path
still receives the reported flexible acceleration. The opt-in
`orifice_force_path="published_owc"` uses those published signal routes and
reduces the 130 s heave and flexible
position differences from `1.94` m and `0.280` m to `0.0366` m and `0.00533` m.
This identifies the primary cause of their divergence. Python keeps the piston
reaction in the flexible state equation by default. The source-compatible
setting is for trajectory comparison; it does not restore the omitted piston
force or establish physical validity above the source's Mach threshold. Surge
and unwrapped pitch still differ by up to `0.174` m and `0.0107` rad, so exact
full-trajectory numerical parity remains unestablished.

The targeted [RM3 sea-state matrix export](https://github.com/cmudrc/wec-sim-python/actions/runs/37682612044)
confirms that MATLAB's applied added-mass matrix and adjusted rigid mass sum
to the original rigid mass plus HDF5 infinite-frequency added mass in the
active surge/heave/pitch block. The paired matrix gate passes for both bodies.
The split preserves the effective active mass exactly, so changing static
added-mass coefficients would have targeted the wrong mechanism.
The ordinary RM3 convolution solver retains implicit added mass. Paired
MATLAB tests also exercise `added_mass_scheme="simulink_delay"` for source
numerics; this compatibility setting is not the default dynamics path.
For Cases 3 and 4, paired gates now run the implicit default independently
of the explicit source-numerics option. The default has maximum float-surge
differences of 56.9 and 48.3 mm against MATLAB; the source option reduces
them to 28.0 and 24.8 mm. This comparison does not change the default
hydrodynamic coefficients or added-mass treatment.
The applied added-mass force in all three saved sea states uses acceleration
extrapolated from the two preceding 0.1 s samples to the current time minus
the Simulink block's `1e-7` s Transport Delay. After accounting for the
postprocessed pitch-inertia correction, that reconstruction matches the
reported force within `7e-7` N or N m across both bodies. This establishes
the source's acceleration feedback at output times. The optional
`simulink_delay` scheme applies the same mass split and delayed feedback.
The exact MATLAB joint adds `PTO stroke * sin(pitch)` to the two bodies'
relative surge, a term absent from the reduced Python geometry. In the saved
MATLAB traces this term peaks at 0.304, 0.222, and 0.140 m across the three
sea states, compared with 0.039 m in regular MCR case 8. Exact geometry alone
regressed regular-wave gates. Combining it with the source's delayed mass
feedback passes both the imported sea states and the previously paired
regular-wave cases. The delay is a Simulink numerical setting, not a measured
WEC property, and is confined to explicitly selected RM3 convolution
comparisons here.

The published RM3 body-to-body Cases 5 and 6 use the same suspect fitted
state-space radiation as the fourth `Radiation_Force_Options` setting. At 400 s,
MATLAB float surge is +0.515 m in Case 5 versus +0.110 m in convolution Case 3,
and +0.686 m in Case 6 versus +0.075 m in convolution Case 4. The fit's
effective common-surge damping is negative at low frequency although source
BEM damping is nonnegative at its sampled frequencies. Matching the MATLAB
state-space trajectory is a diagnostic reproduction, not evidence of physical
validity; it is excluded from validated parity until the fit is resolved.
The [paired Cases 5/6 force run](https://github.com/cmudrc/wec-sim-python/actions/runs/37936907603)
now exports all six MATLAB radiation-force channels for both bodies over the
published 400 s, 0.1 s cases. Given only their saved body velocities, Python
replays the pinned HDF5 state-space fit with cross-body radiation off in Case 5
and on in Case 6. Against those source forces, the largest error is 0.0554% of
a source-channel peak (0.0526% in Case 6); the paired gate is 0.07% with a
3 N peak floor. This tests the force law on prescribed source motion, not an
independent coupled trajectory or the physical validity of the fit. The
default Python convolution dynamics remain unchanged.
The published Cases 5 and 6 change only the state-space radiation flag from
their respective convolution cases, apart from formatting of the inputs.
Over an 8 s window ending at 392 s, MATLAB float-surge means are 0.465 m
and 0.597 m for Cases 5 and 6, compared with 0.101 m and 0.047 m for
Cases 3 and 4. Repeating Cases 5 and 6 with a 0.05 s step gives means of
0.451 m and 0.592 m, so halving the published step does not remove the drift.
At zero frequency the fitted common-surge damping applied by MATLAB is
−14.1 kN s/m without cross-body coupling and −17.5 kN s/m with it.
The pinned MATLAB `bodyClass` sets the state-space direct term to zero even
though the BEMIO HDF5 stores nonzero terms. Including those saved terms in a
diagnostic transfer calculation reduces the negative zero-frequency values
to −3.01 and −3.45 kN s/m, but does not make the fits passive. This is a
low-frequency fit problem consistent with the drift, not proof that the
full nonlinear trajectory has only one cause.
The saved MATLAB radiation-force and velocity histories also constrain that
interpretation: integrating the logged resisting-force dot velocity over all
12 body coordinates gives **positive** 71.7 MJ in Case 5 and 70.7 MJ in
Case 6 (8.3 and 10.7 MJ from surge alone). Thus radiation removes net energy
over each published oscillatory trajectory even as the surge offset grows.
Negative zero-frequency damping identifies a problematic slow mode; it does
not mean the fitted radiation supplies net energy over these 400 s runs.
The diagnostic fixture and a freshly BEMIO-generated HDF5 from the pinned
RM3 Applications input have identical 260 frequency samples and source BEM
damping. Their common-surge fitted transfer curves differ by at most
`6.9e-9` N s/m and projected active-mode minimum damping by at most
`3.0e-9` N s/m. The raw fitted A/B/C arrays differ because they encode
equivalent state-space realizations; the paired gate compares the transfer
responses. Thus the negative low-frequency fit is present in the current
MATLAB case input as well as the historical fixture.
The diagnostic also projects the BEM and fitted radiation matrices onto
the RM3 joint's four moving coordinates: common surge, float heave, spar
heave, and shared pitch. Pitch is expressed as travel at a 20 m lever for
an energy-conjugate matrix with consistent units. The fitted matrix's
least-damped direction is −14.3 kN s/m at zero frequency with cross-body
radiation off and −17.6 kN s/m with it on. The source BEM matrix has small
negative eigenvalues too (minimum −0.422 kN s/m across its sampled
frequencies), so these data do not establish global BEM passivity. At the
sampled frequency nearest the 8 s wave, the coupled source minimum is
−0.030 kN s/m while the fit is −5.87 kN s/m. These are properties of the
linearized radiation data in the active joint subspace, not proof that every
negative mode is excited in Cases 5 and 6.
Diagnostic PR #6 was closed without merging its state-space solver. The
convolution results for Cases 3 and 4 were identical, sample for sample,
between the original branch point and that diagnostic branch.
The pitched-slider map used by the current default solver is also checked
independently of MATLAB output: finite differences of body position verify
its velocity Jacobian and acceleration bias, and an unforced two-body joint
with symmetric added mass conserves kinetic energy within `1e-8` relative
over 10 s. These checks guard the default joint mechanics; they do not
establish passivity of the published fitted radiation model.
An independent no-wave check also separates the models: from a 0.01 m/s
common-surge perturbation, the default convolution solver ends at 0.00959 m/s
without cross-body radiation and 0.00934 m/s with it after 400 s. The pinned
fit's zero-wave linearization instead grows to 0.077 and 0.115 m/s, with
positive real growth rates of 0.0053 and 0.0064 s⁻¹. This is a regression
check on the default dynamics, not a Case 5–6 parity claim. The explicit
`simulink_delay` option is used only for source-numerics comparisons and does
not enable state-space radiation.

The paired case-runner checks compare the MATLAB and Python time grids, all
active body positions and velocities, and stationary degrees of freedom. RM3
checks also compare PTO stroke,
speed, force, and mechanical power; OSWEC checks compare PTO angle, angular
speed, torque, and mechanical power. MATLAB reports absorbed PTO power with a
negative sign, while the Python API reports positive absorbed power. The
Sphere checks now require a stationary 0 m case and narrower motion and force
limits for the displaced cases. Tolerances remain explicit in
`tests/test_case_dynamics_parity.py` and are calibrated against the pinned
MATLAB reference output, not against historical Python fixtures.

Source comparisons: [MATLAB wave class](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/waveClass.m),
[MATLAB wave-number function](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/functions/BEMIO/calcWaveNumber.m),
and [MATLAB body class](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/bodyClass.m).

## Reference model and case inventory

`tests/reference_cases.csv` inventories 63 published input files at the
pinned revisions: two core examples and 61 WEC-Sim Applications cases. It
identifies 21 RM3, 13 OSWEC, and 12 Sphere cases by explicit hydrodynamic
file paths, plus 17 other or dynamic cases. Generate the CSV with
`python tools/build_reference_case_inventory.py CORE_CHECKOUT APPLICATIONS_CHECKOUT`.
This is a source inventory, not a claim that every case runs in Python.

The [latest full reference-model sweep](https://github.com/cmudrc/wec-sim-python/actions/runs/38108295073)
completed 40 MATLAB/Python jobs successfully on `2d35cfb5`, with one optional
comparison skipped. It includes the current-profile comparison. The three-seed
continuous-yaw control was checked in the separate focused CI runs linked
above; the RM3 Cases 5/6 fitted-force audit also uses its separate focused
job. The sweep covers the cases selected by its matrix, not all 63 inventory
entries or the named gaps below.

### Published-case limits and open gaps

This table records both open published-case comparisons and the limits of
cases with passing independent Python motion gates. `Desalination`, WaveStar
Fault Implementation, `FloatingOWC`, and `MOST` have passing coupled-motion
gates for their pinned published configurations; their rows identify outputs
or other regimes that remain unpaired. The full native motion gates still
open here are RM3 body-to-body Cases 5/6 and the state-space radiation option,
OSWEC PassiveYawRegression and the published variable-hydrodynamics irregular
yaw run, RM3 MoorDyn ParaView, and WaveStar NMPC. The ode45 nonlinear-hydro
and OWC OrificeModel physical-force paths have source-timing or force-routing
limits described below. Numerical `Wave_Markers` elevation is paired, while
glyph rendering is not. A passing upstream MATLAB application test only
establishes that WEC-Sim ran; it does not establish Python parity.

| Published case(s) | Evidence and remaining work |
| --- | --- |
| RM3 `B2B_Case5` and `B2B_Case6`; state-space setting in `Radiation_Force_Options` | MATLAB trajectories, fitted-transfer diagnostics, and source-motion radiation-force replays exist for all three settings. The pinned fit has negative low-frequency damping in the active joint coordinates. No physically validated Python state-space trajectory pair exists. |
| OSWEC `PassiveYawRegression` | Three phase realizations of the published irregular hold show up to `83.6%` disagreement in 250 s PTO work after event divergence, despite source-force replay matching all six channels. Prescribing only the source heading-update sequence restores tight yaw and PTO-work gates across all three seeds. In a pinned three-step audit, MATLAB PTO work changes by `-47.8%` then `+31.9%` for seed 1 and `-19.2%` then `+102.0%` for seed 3 as `dt` halves from `0.01` to `0.005` to `0.0025` s, with identical incident waves. Native full-trajectory parity and a converged held-heading reference remain unestablished. |
| OSWEC `Variable_Hydro/Passive_Yaw`, published 600 s irregular sea | The [R2025b published run](https://github.com/cmudrc/wec-sim-python/actions/runs/38035515976) and Python agree on the sea, saved-state heading selection, and six excitation components, but native heading events diverge after 54.73 s. The 600 s native Python run differs by up to `0.1103` rad yaw and `9.05%` PTO work. Prescribing MATLAB's saved heading choices reduces those to `0.00868` rad and `0.0081%`; that is a diagnostic, not native trajectory parity. |
| RM3 `Paraview_Visualization/RM3_MoorDyn_Viz` | The pinned 80 s R2025b application now exports its Traditional 1000-bin sea, both bodies, PTO, MoorDyn connection, and three fairlead tensions in a [focused source comparison](https://github.com/cmudrc/wec-sim-python/actions/runs/38084526959). On multiple independent phase realizations, Python wave elevation agrees within `1.8e-14` m and all body excitation components within `3.7e-8` N or N m. A derived 10 s run supplies every source MoorDyn connection pose and speed; replaying those through the same pinned native library reproduces all six source load components within `0.01` N or N m on both Linux and macOS. The same seeded MATLAB source motion and mooring forces are identical at all shared samples with `dtOut=0.01` and the published `0.1` s; the output interval does not explain the gap. Published-step full coupled motion remains unpaired: the original 80 s case has shown up to `234` mm float surge error across repeated random seas, while independent seed-1 and seed-2 derived seas show `191` and `20.1` mm over 10 s despite matching waves and native mooring loads on saved poses. Halving Python's physical time step from `0.01` to `0.005` s changes seed-1 surge by at most `0.35` mm, far below its source gap. A [derived full-duration fine-step gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38108596566) now pairs 80 s active body and MoorDyn motion, connection loads, and PTO force: maximum float-surge and MoorDyn pitch-moment differences are `0.827` mm and `1.636` kN m against `2` mm and `2` kN m gates. The original published-step gate remains red, and the physical Python default is unchanged. Three actual published wave VTP frames at 0, 10, and 80 s match Python in a [focused R2025b gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38105960778). Six actual float/spar body VTP frames at those times also pass a [pinned prescribed-pose gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38111446848): every facet is matched one to one, with maximum center difference `1e-5` m and area difference `3.50e-5` m². These body frames use MATLAB-saved poses. The [full 801-frame MoorDyn VTP gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38112549662) writes all three lines from saved source node and segment-tension records: connectivity and offsets match exactly, node coordinates within `1e-5` m, and tensions within `0.05` N across all 2,403 line pieces. This is prescribed-source-line output parity; the full body and wave VTP histories and original-step coupled motion remain unpaired. |
| OSWEC `Desalination` | The pinned 300 s Simscape Fluids application runs under R2025b; [expanded source baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37838787789) records flap, PTO, pressure, flow, and mechanical-power traces. Its 250-bin incident sea, six flap excitation components, and body-local rod motion pair with Python over all 30,001 samples; maximum excitation error is below `1e-6` N or N m. The source Morison drag law matches 605 sampled flap states within `1e-6` N or N m after accounting for the source's logged-force sign. The source PTO actuation force exactly follows the measured cylinder force with a one-step delay over the full trajectory. Given MATLAB's inlet pressure, the Python osmotic valve and resistance reproduce all 30,001 permeate-flow values within `2.1e-12` m³/s. Given the source accumulator flow and startup pressure, the Python gas-volume and hard-stop law reconstructs all 30,001 pressure values within `1e-3` Pa using the source's backward Euler flow step. The initial inferred liquid volume is negative because the source starts its network far below the 3 MPa precharge while allowing finite hard-stop penetration. Given the source chamber pressures and rod speed, the Python ideal cylinder reproduces all 30,001 rod-force samples within `1e-7` N and both port-flow traces within `1e-12` m³/s; the measured source force has the opposite sign. The published hydraulic junction balance closes within `2e-15` m³/s, including the relief branch. The high-pressure network now advances from saved rod speed alone, without measured source pressure or branch flows. Over all 30,001 samples it gates pressure within `500` Pa (source pressure reaches 5.6 MPa), permeate, brine, and recovered flow within `1e-5` m³/s each, and relief/accumulator flow within `3e-3` m³/s each. Using the source feed flow instead, the network reconstructs the first 6001 pressure states within `0.01` Pa. This also provides a stronger prescribed-rod-motion component gate than the coupled check below. A four-valve cylinder model using the published passive-orifice settings predicts chamber pressures and rod force from Python high pressure and the saved rod speed. The pinned legacy solver alternates raw chamber pressure at the 0.01 s sample rate while valve flow remains smooth; over 20% of valve 1 samples log flow against the pressure drop. Pointwise chamber pressure and force are therefore not claimed as paired. After a 0.1 s mean, the independent rod-force error is below 100 kN RMS; net rod work differs by less than 3% over 300 s. Raw rod-force disagreement remains above 1 MN RMS. The independent 300 s Python runner now couples incident sea, flap motion, radiation memory, source-convention Morison drag, four-valve chamber force, and high-pressure dynamics. With the pinned HDF5 and saved MATLAB random phases, but no source forces or state trajectories as inputs, its paired gates are maximum flap-pitch error below `0.02` rad, pitch-speed error below `0.01` rad/s, body-center error below `0.1` m, rod-speed error below `0.02` m/s, high-pressure error below `40` kPa, rod-stroke error below `0.05` m, 0.1 s mean rod-force error below `150` kN RMS, and net rod-work error below `2%` over 300 s. The source-specific Morison convention and one-step PTO actuation delay are confined to this case runner. The raw chamber pressure/force artifact remains unpaired. The source also permits chamber gauge pressure below −15 MPa, so this numerical comparison is not evidence of physical cavitation behavior. The two published OSWEC PTO-Sim hydraulic crank applications are covered above. |
| OSWEC `Paraview_Visualization/OSWEC_NonLinear_Viz`; RM3 and OSWEC `Wave_Markers` rendering | Python writes regular-wave VTP surfaces with the source 2D grid, quad connectivity, and `ground.txt` metadata; a two-frame, six-vertex pinned MATLAB R2025b gate compares point coordinates and polygons. A second source-writer gate compares a translated and XYZ-rotated triangular body mesh, cell areas, and all three pressure arrays over two frames. The published OSWEC nonlinear visualization case now has a full 120 s pressure pair on its actual 1,042-facet flap STL and 1,201 saved source poses. Python independently regenerates the regular sea and evaluates facetwise elevation, hydrostatic pressure, and nonlinear and linear Froude–Krylov pressure. Maximum MATLAB/Python differences are `4.4e-11`, `1.8e-10`, and `1.7e-10` Pa, respectively, against a `1e-8` Pa paired gate. The source pose is an input to this pressure check; the independently advanced nonlinear flap motion is paired against a refined MATLAB run below. The full published VTP scene remains unpaired. The 25 numerical elevation histories in each published Wave Markers case are paired above; their rendered marker glyphs remain unpaired. The actual RM3 MoorDyn line VTP history has the prescribed-source-line gate above. |
| WaveStar `WECCCOMP/WECCCOMP_Fault_Implementation` | [Fresh R2025b source runs](https://github.com/cmudrc/wec-sim-python/actions/runs/37871185462) complete the pinned 141.2 s application and export five bodies, four PTOs, the 500-bin sea, hydrodynamic forces, PTO actuation, and the controller's true and faulted position at every 0.001 s step. Python reproduces the saved sea and all six excitation channels within `1e-12` in their units. The float and arm retain the base case's one-coordinate linkage: body positions, speeds, and B-to-C PTO stroke/speed match source states within `1e-12`. The public `StribeckFriction` law reconstructs all three joint-friction torques within `1e-12` N m, including the pivot's 55–115 s fault window. Given the realized position-sensor disturbance and source stroke, `WaveStarFaultController` reconstructs axial motor force within `1e-7` N after the published discrete filters and actuator transfer. The measured sensor has approximately 3% exact dropouts and 0.003 rad noise. The public `run_wavestar_fault_published` uses only that exogenous stochastic disturbance and the saved wave phases as inputs; it independently advances its own 92 radiation states, arm angle, friction, sensor geometry, controller state, and PTO force at 0.001 s. Across two distinct random realizations and all 14,121 output samples, its worst float-pitch error is `9.16e-5` rad, pitch-speed error `0.00702` rad/s, float-position error 44 µm, PTO-stroke error 18.3 µm, axial-force error 0.128 N, and radiation-force error 0.0026 N or N m. The paired gates allow up to `2e-4` rad pitch, `0.015` rad/s speed, `0.1` mm body position, `0.04` mm PTO stroke, `0.3` N axial force, and `0.01` N or N m radiation force. These are numerical gates for the source's fitted radiation and sampled controller; the broader `WEC` radiation default remains unchanged. |
| WaveStar `WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive` | A [fresh R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37876181993) completes the published Sea State 6 application with its 225 s `ode8` simulation, autoregressive predictor, and NMPC. The selected source gate pairs its 500-bin JONSWAP sea, six float excitation channels, rigid float/arm linkage, B-to-C PTO stroke and speed, and the 10–15 s resistive startup at full 0.05 s output resolution. Given the logged torque command and rod stroke, public `WaveStarNmpcActuator` independently reproduces axial PTO force within `5.7e-11` N across all 4,501 samples. Its controller kinematic constants differ slightly from the physical rod geometry, so the two moment arms are kept separate. The source predictor briefly reaches `5.24e6` N m during 10.05–10.95 s before NMPC activates at 15 s; this is a source trace, not a Python model result. Given the saved motor stroke and preceding torque command, public `WaveStarNmpcObserver` independently replays its controller-side inverse linkage, sampled velocity filters, and five-state Kalman update within `7e-13` across all 4,501 samples. Given the saved estimated excitation moment, public `WaveStarNmpcPredictor` reproduces the source [forward-backward AR(18) fit](https://www.mathworks.com/help/ident/ref/ar.html) and all 40 forecast steps. Before 10 s both are zero. During the source's unstable 10–11 s startup fit, the maximum absolute forecast difference is `0.056` N m while the source reaches `5.24e6` N m; after the 11 s refit, every remaining value through 225 s agrees within `1e-9` N m. Composing the independent observer and predictor from saved stroke and preceding torque request retains that `1e-9` N m gate after 11 s; in the unstable 10–11 s window, tiny observer roundoff amplifies to at most 25 N m against the source forecast of up to `5.24e6` N m. The public `WaveStarNmpcController` replays the source RTI quadratic program, including the zero-command period and resistive startup. Given the saved estimator states and excitation forecasts, all 4,501 source torque commands agree within `1e-4` N m (maximum `6.25e-5` N m; RMS `1.63e-6` N m). Composing the independent observer, predictor, and controller from saved PTO stroke and preceding torque request also meets the `1e-4` N m gate (maximum `6.21e-5` N m). The MATLAB stroke and preceding command are still inputs to this replay; independent closed-loop motion remains unpaired. Before any PTO actuation, the published 0.05 s `ode8` source differs from fine-step unforced Python motion by at most `0.00342` rad pitch, `0.0354` rad/s pitch speed, `1.67` mm float-center position, `0.681` mm PTO stroke, and `0.180` N or N m radiation force across the 200 samples through 9.95 s. A [derived fine-step MATLAB run](https://github.com/cmudrc/wec-sim-python/actions/runs/37879718043) keeps the same 500-bin sea, geometry, hydrodynamics, and zero PTO actuation, but uses `ode4` at 0.001 s through 9.95 s. Its saved wave components match the published run exactly. Across all 9,951 samples, independent Python motion is gated within `2e-5` rad pitch, `4e-4` rad/s pitch speed, `10` µm float-center position, `0.2` mm/s float-center speed, `4` µm PTO stroke, `0.08` mm/s PTO speed, and `5e-4` N or N m radiation force. Observed maxima are `1.56e-5` rad, `3.58e-4` rad/s, `7.60` µm, `0.175` mm/s, `3.11` µm, `0.0713` mm/s, and `4.51e-4` N or N m. At matching 0.05 s samples, the published coarse MATLAB run differs from the fine MATLAB run by `0.003414` rad pitch and `0.03530` rad/s speed, nearly the entire earlier Python gap. This isolates that pre-control gap to the source numerical path (step size and/or solver). A [derived fine-step resistive-control run](https://github.com/cmudrc/wec-sim-python/actions/runs/37882941804) uses the same pinned sea, geometry, and hydrodynamics with `ode4` at 0.001 s through 14.95 s; its controller and actuator also sample at 0.001 s, and the unused AR predictor is deferred to 15 s. On saved stroke and preceding command, the Python observer reproduces all five source states within `1e-10`, and the actuator reproduces axial force within `1e-8` N. An independent Python closed loop from wave phases and its own motion is gated over all 14,951 samples within `5e-5` rad pitch, `4e-4` rad/s pitch speed, `25` µm float-center position, `0.25` mm/s float-center speed, `10` µm PTO stroke, `0.1` mm/s PTO speed, `0.005` N m command, `0.03` N axial force, and `0.001` N or N m radiation force. Observed maxima are `3.75e-5` rad, `3.58e-4` rad/s, `18.5` µm, `0.175` mm/s, `7.45` µm, `0.0713` mm/s, `0.00362` N m, `0.0206` N, and `0.000788` N or N m. This establishes independent motion with resistive control on a matched fine-step path. Public `WaveStarNmpcPTO` now composes the observer, AR predictor, RTI controller, and actuator with a 0.05 s control step and a 0.001 s plant step. An independent 225 s Python run using only the pinned sea phases and hydro data completes all 4,501 controller samples with finite bounded motion and commands. Against the published 0.05 s `ode8` source, its maximum differences are `0.0392` rad pitch, `0.414` rad/s pitch speed, `7.82` mm PTO stroke, `2.98` N m torque command, and `14.87` N axial force. Pitch already differs by `0.00342` rad before control starts, and the gap grows during feedback. These full-duration values are a diagnostic, **not a parity gate**; the published closed-loop NMPC trajectory remains unpaired. |
| OWC `FloatingOWC` | The pinned [source and air-train run](https://github.com/cmudrc/wec-sim-python/actions/runs/37895804527) completes the published 500 s `ode45` case and saves 50,001 samples of both bodies, PTO, wave, six-component MoorDyn coupling, available fairlead tensions, turbine speed/control, efficiency, and pneumatic/turbine power. BEMIO generates its two-body, 12-DOF HDF5. Python now pairs the regular-wave excitation on all 12 body channels over all 50,001 samples; the largest absolute difference is `4.04e-6` N or N m, below a `1e-5` gate. This checks hydrodynamic input only, with no source force history passed to Python. The published run leaves body-to-body coupling disabled, so each body applies `6×6` added-mass and damping blocks despite the shared HDF5; Python matches those logged matrices within `1e-5` and `1e-7` respectively. A corrected HDF5 hydrostatic orientation matches both full `6×6` source matrices within `1e-5`, including asymmetric yaw couplings. The source radiation-force record also equals its damping matrix times body velocity within `1e-8` N or N m. The published nine-line configuration logs five named `FairTen` channels. The state-bus exporter records chamber pressure, water-column displacement and speed, and rotor speed; Simulink names the bus speed field `signal4`. Public `FloatingOwcColumnJoint` now reconstructs the water-column center pose and velocity from the saved floater state and PTO axial stroke/speed at all 50,001 samples within `1e-10` m, m/s, or rad/s. Its world heave relative to the column equilibrium center and world heave speed reproduce both chamber input signals within `1e-10`. This is a prescribed-floater-and-PTO kinematic gate. A separate two-coordinate Python heave solve advances both bodies from the pinned wave and HDF5 coefficients using only the saved MoorDyn vertical force and chamber pressure as external histories: over 500 s, maximum floater/column heave differences are `0.00516`/`0.00981` m and heave-speed differences are `0.00258`/`0.00355` m/s (gates `0.015` m and `0.006` m/s). This reduced model omits surge and rotation and prescribes rather than advances MoorDyn and chamber state; it is a mechanical-motion diagnostic, not full floating OWC parity. Given only the saved water-column displacement and speed after initialization, public `FloatingOwcChamber` and `FloatingOwcTurbine` advance the coupled chamber pressure and rotor speed over all 50,001 samples. Local ARM64 maximum differences are `0.0121` Pa pressure and `0.000501` rad/s speed; RMS differences are `0.00416` Pa and `0.000263` rad/s. Linux x86 gives `0.0124` Pa pressure RMS on identical source data. The explicit gates are `0.05` Pa and `0.002` rad/s maximum, `0.02` Pa and `0.001` rad/s RMS. From locally advanced states, maximum differences are `0.0195` W for logged turbine load power, `0.141` W pneumatic power, `5.01e-5` N m load torque, and `9.81e-6` turbine efficiency; power gates allow `0.1` W and `0.5` W. Evaluated directly on saved pressure and speed, the fitted Wells curves and control law reproduce each logged turbine output within `1e-7` in its units. This is a prescribed-water-column-motion air-train comparison. A separate fully coupled seven-coordinate Python trajectory with live MoorDyn, chamber, and turbine feedback is now gated below. |
| Other/dynamic: `MOST` | A focused pinned MATLAB reader gate compares Python's decoding of the checked-in 20,750-step TurbSim `.bts` file at selected grid points and with spatial checks at every time. A second source gate executes the published `RunTurbsim.m` on a 1,000-frame excerpt and compares all 249 output times, three velocity components, ten X planes, and 12×12 YZ grid points with Python's lazy advection view; the maximum difference is exactly zero. A third gate compares the six-component force and all three horizontal/vertical fairlead tensions from the pinned direct `nonLinearStaticMooring.m` at ten translated and rotated poses. Observed maximum differences are below `1e-5` N force, `1.5e-4` N m moment, and `5e-6` N tension; the paired gates are `1e-3` N, `1e-2` N m, and `1e-3` N respectively. A [fourth pinned MATLAB gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37957379862) extracts the active BEM rotor function from `MOST_Lib.slx`, generates its IEA 15 MW blade inputs, and compares Python's three-blade root forces and moments at five hub/wind states, including moving and rotated hubs and spatially varying wind. The maximum observed component difference is `1.5e-8` N or N m, below the absolute `1e-4` gate. The published input sets `nonlinearStaticData.flag=1`, leaves `lookupTableFlag=0`, and sets `aeroLoadsType=1` for BEM; its active mooring block is `MooringNLStatic`. The generated mooring and aerodynamic lookup tables are not used by this case. These gates establish wind input, frozen-turbulence time shift, direct static catenary force, and isolated BEM blade loads. The 10 s rotor gate below advances turbine state from prescribed platform motion. A six-coordinate runner couples it to the platform over two pinned 10 s JONSWAP seas (4 m/seed 1 and 6 m/seed 2); the 30 and 60 s developed-sea gates below extend the 6 m condition beyond its wave ramp; a separate 12 m/s constant-wind gate exercises nonzero initial pitch. The published 1,000 s coupled case now has a local independent trajectory comparison and a [passing full-duration source gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38071234135). Other wind/controller regimes remain unpaired. The WaveStar base and fault applications and both WaveBot applications are paired separately above. |
| Other/dynamic: published-step `Nonlinear_Hydro/ode45` and `OWC/OrificeModel` | Both nonlinear-hydro ode45 cases pass full-duration derived fine-step gates, but the original 0.05 s source trajectories retain force-timing and step-sensitivity gaps. The OrificeModel physical coupled path diverges late, while an explicit source-convention force route passes 130 s body/PTO gates; this does not validate the source air model beyond its Mach threshold. |

The [RM3 MoorDyn visualization force export](https://github.com/cmudrc/wec-sim-python/actions/runs/38085085794)
records acceleration, excitation, radiation, added mass, restoring, Morison,
linear damping, and total hydrodynamic force at every 0.01 s step of the
derived 10 s sea. The signed components reconstruct the reported total to
less than `1e-7` N or N m for both bodies. The source's added-mass output is
the applied force after its Simscape mass adjustment, not simply the original
HDF5 `A_inf` times reported acceleration: at `0.01 s`, the float's reported
surge component is `+2.257 MN`, while original `A_inf * acceleration` is
`-0.0793 MN`. The source coefficient shifts `2.837 million kg` into adjusted
body mass and uses delayed acceleration feedback. Comparing the logged
hydrodynamic-force column directly with Python's physical added-mass term
would therefore misidentify a reporting convention as a dynamics defect.
After undoing the source's rotational-force postprocessing, the adjusted
rigid inertia projected onto the RM3 joint agrees with the logged
hydrodynamic, PTO, and MoorDyn generalized forces across all 1,001 samples:
maximum residual is `1.37e-6` N or N m. At startup the source's shifted mass
and zero delayed-feedback state give float/spar surge accelerations of
`-0.183/-0.00725` m/s². With the original HDF5 added mass implicit, those
same initial applied forces predict `-0.526/-0.200` m/s², close to Python's
first-step mean `-0.523/-0.199` m/s². This establishes the initial source
mass-feedback transient; it does not establish that the transient causes the
entire trajectory gap. On the same 10 s sea, the opt-in output-step
`simulink_delay` approximation has `0.198` m float-surge error versus `0.191`
m for the physical implicit solver. The coupled RM3 trajectory gate remains
red; matching the source's variable-step delayed feedback requires further
evidence before changing Python's physical default.
The reproducible
[source-driven counterfactual](tools/analyze_rm3_moordyn_feedback.py)
projects the difference between
the logged delayed added-mass force and the force from current acceleration
through the four-coordinate joint, then applies that prescribed residual to
Python's physical solver with live MoorDyn. On this seeded 10 s sea, its
float-surge error falls from `0.191` m to `0.136` m when the correction stops
at `0.1` s, `0.0215` m at `0.5` s, and `0.00637` m at `1` s; applying it for
all 10 s gives `0.00407` m. Thus the first second of the source's numerical
mass feedback explains most of this particular motion difference. The
correction depends on the saved MATLAB acceleration and force history; it is
an attribution experiment, **not** an independent Python parity result or a
candidate default force model.
The [second-seed R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/38087382937)
independently reproduces all 1,000 MATLAB random phases from `phaseSeed=2`,
then pairs wave elevation, six-component excitation, source force balance,
and adjusted-mass joint balance. Its 10 s float-surge difference is `20.1` mm
with the physical Python solver. Applying the same source-driven correction
for only the first `0.1`, `0.5`, or `1` s reduces the final difference to
`11.9`, `3.02`, or `1.57` mm; correction for all 10 s gives `1.54` mm.
The separate seed supports startup mass-feedback attribution without
providing an independent replacement for MATLAB's delayed feedback or
closing the published coupled-motion gate.

A [source maximum-step audit](https://github.com/cmudrc/wec-sim-python/actions/runs/38088631075)
holds the seed-1 sea, `ode45`, `dtOut=0.01` s, hydrodynamics, PTO, and native
MoorDyn input fixed while changing only MATLAB `simu.dt`, which sets
Simulink's `MaxStep`. At `0.01`, `0.005`, and `0.0025` s, all 1,000 random
phases and all 1,001 saved wave samples are identical. MATLAB's 10 s float
surge is `0.860672`, `1.083537`, and `1.084102` m; spar surge is
`1.348353`, `1.548801`, and `1.549752` m. Thus the published-step to half-step
change is `223` mm float and `200` mm spar, while the next halving changes
them by only `0.565` and `0.951` mm. Python's physical implicit-mass path at
`0.01` s has float surge `1.051404` m: its difference from the `0.0025` s
MATLAB run is `32.7` mm, versus `190.7` mm from the original `0.01` s MATLAB
run. The physical Python `0.01` to `0.005` s step change was below `0.35` mm
on this sea. The coarse MATLAB source trajectory is therefore strongly
maximum-step-sensitive. These two refinements do not prove convergence of
all output channels or isolate which internal wave, solver, or MoorDyn
substep causes the sensitivity. They also do not explain the remaining
`32.7` mm; the independent
coupled-motion gate stays red, and Python's physical default is unchanged.
The [refined-force gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38089549994)
also saves all six force components on both finer source runs. Their logged
force sums and adjusted-mass four-coordinate balances close to within
`1.2e-6` N or N m, and Python independently regenerates the source phases,
wave, and excitation. Against the `0.0025` s MATLAB trajectory, the physical
Python path differs by `32.7` mm float and `30.4` mm spar surge at 10 s.
Applying the **source-derived** delayed-mass residual to Python for only the
first `0.1` s reduces those final differences to `3.61` and `3.41` mm;
applying it throughout 10 s gives `4.22` and `3.97` mm. The refined case's
remaining difference is therefore also dominated by the source's initial
mass-feedback transient. This correction uses MATLAB's saved acceleration
and force, so it remains an attribution diagnostic, not an independently
generated parity trajectory or a reason to change the physical default.

The independent physical Python solver was also refined on this same
`0.0025` s source sea, using the pinned native MoorDyn library and HDF5 but
no saved source force or motion as input. Halving Python's step from `0.0025`
to `0.00125` s changes either body's position by at most `0.109` mm and speed
by at most `0.032` mm/s over 10 s; a
[paired Linux run](https://github.com/cmudrc/wec-sim-python/actions/runs/38097663337)
confirms the same maxima within rounding. The focused test gates those changes
at `0.25` mm and `0.1` mm/s. The float's final surge difference from the
finer MATLAB source changes only from `32.17` to `32.07` mm. Selecting the
opt-in output-step delayed-mass scheme instead gives `30.29` and `31.13` mm
at the two Python steps, so refining that approximation does not close the
source gap either. These independent runs rule out coarse Python stepping
as the main explanation for the refined-source mismatch; the strict
published coupled-motion gate remains open.

A [third source step refinement](https://github.com/cmudrc/wec-sim-python/actions/runs/38098351547)
on the same 10 s seed-1 sea halves MATLAB's ode45 maximum step again to
`0.00125` s. All saved phases and wave samples still match, and the source
force sum and adjusted-mass joint balance gates pass. The float/spar surge
changes by `37.59`/`34.75` mm relative to the `0.0025` s source, after the
preceding `0.005` to `0.0025` s halving changed them by only `0.72`/`0.95` mm
at most. This nonmonotonic change rules out treating the earlier small
step difference as convergence. Python's independent physical `0.00125` s
run is within `5.52`/`4.85` mm float/spar surge of this latest MATLAB run at
10 s. The source MoorDyn pitch moment changes by up to `37.60` kN m between
the last two MATLAB steps; Python's physical `0.0025` to `0.00125` s change
is at most `1.00` kN m. Against the latest source, the independent Python
moment still differs by `8.79` kN m, above the `2` kN m published-case gate.
The closer motion agreement is not a validated converged pair. The published
80 s coupled-motion gate still fails, and the Python physical default remains
unchanged.

A [fourth source coupling-step refinement](https://github.com/cmudrc/wec-sim-python/actions/runs/38102253800)
halves `simu.dt` to `0.000625` s, changing both ode45 MaxStep and MoorDyn's
call interval. The source sea and logged force balances still agree. Relative
to `0.00125` s, source float/spar surge changes by up to `5.99`/`5.21` mm
and MoorDyn pitch moment by `9.69` kN m over 10 s; source convergence remains
unproven. An independent Python physical run at `0.00125` s, without saved
source motion or forces as inputs, matches this refined source within `0.466`/
`0.360` mm float/spar surge, `5.19e-6` rad pitch, `132`/`104` N MoorDyn
surge/heave force, and `0.906` kN m MoorDyn pitch moment on local ARM64.
Explicit paired gates check both bodies' active positions and velocities and
the MoorDyn connection pose, velocity, and force: `2` mm or `5e-5` rad
positions, `0.5` mm/s or `5e-5` rad/s speeds, and `500` N or `2` kN m loads.
This is a derived 10 s
numerical pair, not convergence proof or full published-case parity; the
physical default and the red 80 s gate are unchanged.

A [fifth and final short-case source refinement](https://github.com/cmudrc/wec-sim-python/actions/runs/38103660538)
halves `simu.dt` to `0.0003125` s on the same sea. The `0.000625` and
`0.0003125` s MATLAB trajectories differ by at most `0.195`/`0.218` mm
float/spar surge, `1.09e-6` rad pitch, and `1.391` kN m MoorDyn pitch
moment. All active body and mooring pose, speed, and load changes pass the
same explicit limits as the Python pair. Against the finest source, an
independent physical Python run at `0.00125` s differs by at most `0.271`/
`0.142` mm float/spar surge, `6.16e-6` rad pitch, `110`/`90` N mooring
surge/heave force, and `0.813` kN m mooring pitch moment on local ARM64.
The gates cover the complete 10 s active trajectories, not only final values.
This supports the derived fine-step physical dynamics over that sea at the
declared numerical resolution. It does not validate every source output or
the original `simu.dt=0.01` s, 80 s published trajectory, whose paired gates
remain unmet.
The published 80 s test keeps its original numerical limits and runs as a
strict expected failure for its known coupled-motion and MoorDyn channels.
Wave agreement and any newly failing channel remain required checks; an
unexpected full pass also fails CI until this limitation is reviewed and the
expected-failure marker is removed. This scopes the passing short-case result
without claiming published-case motion parity.

A [deterministic full-duration comparison](https://github.com/cmudrc/wec-sim-python/actions/runs/38099166656)
runs the published 80 s input with only `waves.phaseSeed=1` added in a
temporary copy. Its phases, wave, both-body motion, and MoorDyn connection
loads agree exactly with the first 10 s of the same-seed short source run
at common samples; the source duration and ParaView setting do not explain
that initial gap. Against this repeatable 80 s source, independent physical
Python reaches `203` mm maximum float surge error and `101` kN m maximum
MoorDyn pitch-moment error, above the existing `10` mm and `2` kN m gates.
The [focused CI](https://github.com/cmudrc/wec-sim-python/actions/runs/38101364860)
has 22 passing checks and this one failing coupled-motion
gate. Seeding makes the failure reproducible but does not establish source
step convergence or full coupled parity.

A [full-duration fine-step comparison](https://github.com/cmudrc/wec-sim-python/actions/runs/38108596566)
keeps that seed-1 published sea and all WEC settings, while turning off
ParaView recording, saving at 0.01 s, and refining `simu.dt` (both ode45
MaxStep and the MoorDyn coupling interval) to 0.000625 and 0.0003125 s.
All 1,000 phases and wave elevation at common samples match the published
run; maximum elevation difference is `1.3e-14` m. Over all 8,001 output
samples, the two refined MATLAB runs differ by at most `0.225` mm in spar
surge, `1.23e-6` rad in pitch, `1.391` kN m in MoorDyn connection pitch
moment, and `19.7` N in PTO force. An independent Python run at 0.00125 s
uses the source phase values but no source motion or force histories. Against
the finest MATLAB run, its maximum float/spar surge differences are
`0.827`/`0.762` mm, pitch difference `7.90e-6` rad, MoorDyn connection
surge/heave-force differences `214`/`181` N, pitch-moment difference
`1.636` kN m, and PTO-force difference `178` N. All active body and
connection positions and velocities also pass the existing fine-step gates:
`2` mm or `5e-5` rad position, `0.5` mm/s or `5e-5` rad/s speed,
`500` N or `2` kN m connection load, and `2` kN PTO force. A local ARM64
run reproduces the Linux errors to the shown precision.

The original 0.01 s MATLAB trajectory differs from the finest MATLAB
trajectory by up to `203.4` mm float surge and `104.2` kN m MoorDyn pitch
moment; Python differs from that original by `203.2` mm and `104.0` kN m.
The common sea and the refined source pair support a numerical-step cause
for nearly the entire published-step gap on this seed. Because refining
`simu.dt` changes both ode45 MaxStep and MoorDyn coupling, it does not
isolate either contribution or prove a unique converged limit. This is
full-duration **derived fine-step** parity for the named channels, not
parity with the original published-step trajectory or all fairlead and
ParaView outputs. The physical Python default and strict expected failure
for the published-step gate remain unchanged.

The same pinned run exports actual published ParaView wave surfaces at 0, 10,
and 80 s. Python independently rebuilds their 2,000-point grids from the
exported sea components. All point coordinates match exactly at the source
VTP precision; polygon connectivity, offsets, and `ground.txt` also match.
This establishes selected published wave visualization output. Selected
body frames and the full prescribed-source-line MoorDyn VTP history have
separate gates above; full wave/body VTP histories and original-step coupled
motion remain unpaired.

For the published OSWEC nonlinear visualization case, the pinned source also
exports flap velocity and all six-component force channels. On all 1,201 saved
MATLAB poses, Python integrates its independently calculated facet pressures
over the 1,042-facet STL. Hydrostatic restoring differs by at most `2e-8`
N or N m, and Froude–Krylov-corrected wave excitation by at most `1e-6`
N or N m. The regular-wave HDF5 radiation damping differs by at most `1e-8`
N or N m. These are explicit paired force gates on prescribed poses. The
source's applied added-mass force, after its reported rigid-inertia term is
removed, follows the two preceding acceleration samples through its
`1e-7` s Transport Delay within `0.01` N or N m. This source timing is
separate from Python's default implicit added mass.

A [four-step R2025b refinement](https://github.com/cmudrc/wec-sim-python/actions/runs/38093454463)
keeps the published 120 s, 2.5 m/8 s regular wave, 40 s ramp, hinge, and mesh
hydrodynamics. Disabling pressure and VTP output changes the original 0.1 s
pitch trace by exactly zero. Halving the MATLAB `ode4` and nonlinear-force
sample steps from `0.1` to `0.05`, `0.05` to `0.025`, and `0.025` to
`0.0125` s changes pitch by at most `0.005581`, `0.001295`, and
`0.000317` rad, respectively. The public `WEC.fixed_hinge` path now
advances its own pressure-loaded flap from the STL, HDF5, wave, mass, and
hinge inputs. It uses implicit added mass and resolves pressure between
its 0.1 s output samples. Against the finest source run, maximum 120 s
errors are `0.000135` rad pitch, `0.000110` rad/s pitch speed, `0.622` mm
surge, `0.262` mm heave, `0.522` mm/s surge speed, and `0.240` mm/s heave
speed. The paired gates are `0.0002` rad, `0.00015` rad/s, `1` mm,
`0.4` mm, `0.8` mm/s, and `0.35` mm/s in that order. Reducing Python's
output step from 0.1 to 0.05 s changes pitch by at most `1.1e-6` rad.
This establishes independent motion parity for this fixed-hinge regular-wave
case against a refined MATLAB path. The published coarse-step motion and
other nonlinear six-DOF layouts have separate numerical and feature limits.

The [pinned MoorDyn ParaView writer gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38080388695)
interpolates two lines with three and four nodes to two intermediate frames.
Python's public `ParaviewClass.write_paraview_vtp_mooring` matches all source
point coordinates within `1e-5` m, segment tensions within `1e-6` N, and line
connectivity and offsets exactly. The source `%i` formatter preserves
fractional tensions in scientific notation; Python also retains those values.
The actual 801-frame `RM3_MoorDyn_Viz` line history is paired above using
saved source nodes and segment tensions. Native Python coupled motion at the
original published step remains unpaired.

The MOST [wind-to-BEM source gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37958579712)
feeds three actual time-shifted frames from the pinned TurbSim excerpt into
the active MATLAB rotor function and the Python `MostWindField.sampler` plus
`MostBEM.loads` path. At interior probe points, the largest interpolated wind
component difference is `8.9e-16` m/s (gate `1e-12` m/s). Across all three
blades and six root-load components at each frame, the largest difference is
`7.5e-9` N or N m (gate `1e-4`). This pairs the wind-to-rotor load path for
sampled hub states. The coupled turbine controller and time-evolving platform
are paired separately in the full-duration gate below.

The pinned `MOST_Lib.slx` wind `From Workspace` blocks specify `Holding final
value` after their last input sample. The checked-in TurbSim record yields
19,999 advected frames, with its last input at approximately 999.90 s, while
the published application ends at 1,000 s. `MostWindField.sampler_at` now
holds that final spatial frame for later finite times, matching the source
input convention at the endpoint. The full-duration coupled comparison below
includes that endpoint.

The MOST [baseline-controller source gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37960749134)
regenerates the published IEA 15 MW steady-state table and control parameters,
then runs the active Simulink Baseline block with a prescribed rotor-speed
history for 40 s at 0.01 s. Public `MostBaselineController.iea15mw()` agrees
with all 4,001 source generator-torque and blade-pitch samples within
`1e-4` N m and `1e-9` rad; local maximum differences are `5.6e-8` N m and
`2.3e-15` rad. Both source output traces reach their respective rate limits.
This pairs the isolated controller given rotor speed; rotor and platform states
are not independently advanced by this gate.

A [pinned MOST short coupled run](https://github.com/cmudrc/wec-sim-python/actions/runs/37963619701)
retains the published BEM rotor, direct nonlinear-static mooring, irregular
waves, TurbSim wind, and baseline controller, changing only the end time to
10 s, output interval to 0.01 s, and explorer setting. Given its saved
six-DOF platform position and velocity as the moving support, public
`MostRotor` independently advances rotor speed, azimuth, and controller
state from the checked-in wind excerpt and raw blade tables. Over all 1,001
samples, local maximum differences are `0.000227` rpm rotor speed,
`0.000139` rad azimuth, and `711` N m generator torque; the paired gates
allow `0.001` rpm, `0.001` rad, and `2,000` N m. Pitch stays at zero in this
early, below-rated interval; the separate controller gate above exercises
pitch and torque rate limits. Python wind at the logged hub-height point
differs by at most `8.9e-6` m/s (gate `1.5e-5`). The rotor's own BEM force,
controller state, and azimuth evolve without saved turbine state as input.
Platform motion is still prescribed, so this is one-way rotor trajectory
parity, not a fully coupled MOST platform/turbine trajectory.

The expanded [MOST platform-force source trace](https://github.com/cmudrc/wec-sim-python/actions/runs/37966607366)
records all six hydrodynamic force channels, mooring position and load,
tower-base load, body acceleration, and realized JONSWAP sea over the same
1,001 samples. Python regenerates all 500 frequencies, spectral amplitudes,
bin widths, and phases to `7e-15` or better; wave elevation differs by at
most `2.8e-14` m and six-component excitation by at most `2.3e-6` N or N m
(gate `1e-5`). The direct Python catenary evaluated at every logged moving
mooring pose differs by at most `2.5e-5` N and `4.9e-4` N m (gates `1e-3` N
and `1e-2` N m). The source body force channels add exactly as excitation
minus radiation, added mass, restoring, viscous drag, and linear damping.
The public `MostPlatformHydrodynamics` evaluates the remaining three
state-dependent hydrodynamic terms on the source trajectory from pinned
VolturnUS HDF5, platform mass properties, and the published drag matrix.
Across all 1,001 samples, hydrostatic restoring differs by at most
`3.0e-8` N or N m, quadratic drag by `5.9e-11` N or N m, and the 60 s
radiation convolution by `2.1e-8` N or N m (paired gates `1e-5`, `1e-8`, and
`1e-6`, respectively). These are prescribed-state force checks.
Using the pinned platform mass and infinite-frequency added mass, then
rotating the logged tower-base force from platform to world coordinates,
closes the source's three translational force balances within `1.4e-8` N
(gate `1e-4` N) at every sample. This fixes the force signs, frame, and
Simscape added-mass shift needed for an independent platform solve; the
tower reaction and platform trajectory are still source inputs to this
diagnostic. The separate Newton–Euler reaction gate below resolves the
six-axis tower load; a six-coordinate coupled trajectory is gated below.

Public `MostTowerReaction` now calculates that tower-base wrench from the
published tower, nacelle, yaw bearing, hub, and three blade mass properties,
gravity, rotor and platform kinematics, generator torque, and blade-root
aerodynamic loads. Across all 1,001 states of the pinned 10 s MOST run, its
maximum force-component error is below `1.5e-8` N. Moment errors are below
`0.033` N m about the tilted shaft axes and `1.5e-7` N m in the dominant
pitch moment; the paired gates allow `1e-5` N and `0.1`, `1e-4`, `0.01` N m
respectively. This gate takes MATLAB platform acceleration, rotor state,
generator torque, and BEM root loads as inputs. It verifies the physical
reaction law and its frames. The live turbine/platform trajectory
below removes those saved turbine inputs.

The public `MostPlatformHydrodynamics.simulate` separately advances platform
surge, heave, and pitch from the Python JONSWAP excitation, live catenary,
VolturnUS radiation and added mass, hydrostatic restoring, and quadratic drag.
Only the six-component MATLAB tower-base reaction is prescribed, in the
platform frame; the platform trajectory is independently advanced. Against
the 1,001-sample, 10 s source run, maximum position errors are `0.1435` mm
surge, `0.0716` mm heave, and `6.38e-6` rad pitch; maximum speed errors are
`2.84e-5` m/s surge, `4.42e-5` m/s heave, and `3.93e-6` rad/s pitch.
The paired gates are `0.3` mm, `0.15` mm, `1e-5` rad, `7e-5` m/s, and
`1e-5` rad/s. Sway, roll, and yaw are omitted in this reduced solve, and
the tower reaction is supplied rather than computed from live Python turbine
states. This establishes one-way platform dynamics, not full MOST coupling.

Public `MostCoupled` now starts from the equilibrium platform and independently
advances the Python BEM rotor, baseline controller, six-coordinate platform,
and Newton–Euler tower reaction from only the pinned wind and Python-generated
wave excitation. The tower's acceleration-dependent reaction is assembled as
implicit platform inertia; the whole-trajectory rotor/platform iteration
converges on this 10 s case in nine passes, with the final platform-position
change below `1e-8` m or rad. No MATLAB motion, rotor, blade-load, platform
acceleration, or tower-load history enters the solve. Against the pinned
1,001-sample source run, local maximum position differences are `0.114` mm
surge, `0.287` mm sway, and `0.078` mm heave; maximum roll, pitch, and yaw
differences are `2.28e-5`, `2.73e-6`, and `9.54e-5` rad. Rotor-speed
difference is `0.000164` rpm. The paired gate checks all six positions and
speeds, rotor speed, azimuth, generator torque, and convergence. This
establishes bidirectional six-coordinate platform/turbine motion for the
published 10 s, below-rated case. A second pinned condition changes only the
JONSWAP significant height from 4 to 6 m and phase seed from 1 to 2. The
independent Python solve converges in eight passes; the saved wave bins,
amplitudes, and phases match the fresh MATLAB source to `1e-11`. Wave elevation
and all six excitation-force components differ by at most `2.9e-14` m and
`6.7e-6` N or N m, below paired gates of `1e-12` m and `1e-4` N or N m.
Across its 1,001 samples, maximum platform position differences are `0.114` mm surge,
`0.280` mm sway, `0.078` mm heave, `1.85e-5` rad roll, `2.70e-6` rad pitch,
and `9.05e-5` rad yaw. Maximum velocity differences are `0.0231` mm/s surge,
`0.0952` mm/s sway, `0.0442` mm/s heave, `9.34e-6` rad/s roll, `5.94e-7`
rad/s pitch, and `2.67e-6` rad/s yaw. Rotor speed, azimuth, and generator
torque differ by at most `0.000163` rpm, `1.09e-4` rad, and `502` N m. The
yaw-velocity gate is `3e-6` rad/s across both seas; no dynamics coefficients
were changed for this condition. Both cases cover only the first 10 s of a
20 s wave ramp with the same below-rated wind. The 30 and 60 s developed-sea
extensions below cover motion after the ramp. A separate 12 m/s constant-wind
case below tests the controller with nonzero steady-state initial pitch.
The published 1,000 s trajectory is paired below; other wind/controller
regimes remain unpaired.

The 10 s gates use each MATLAB run's generated controller and steady-state
tables. Two fresh 6 m source runs started 0.0164 rpm and 62.8 kN m apart in
rotor speed and generator torque, and their surge traces separated by
0.787 mm at 10 s. Pairing those generated inputs restores the existing
trajectory gates without changing rotor or platform dynamics.

A [pinned 30 s MOST developed-sea run](https://github.com/cmudrc/wec-sim-python/pull/140)
extends the 6 m, seed-2 condition without changing its wind speed or
controller scripts. Its 1,500-frame TurbSim excerpt produces 749 advected wind frames;
Python reproduces all three components on every saved 10×12×12 frame within
`1e-12` m/s. Over all 3,001 samples, Python's independently synthesized wave
elevation and six-component excitation differ from MATLAB by at most
`1.8e-13` m and `2.5e-5` N or N m. Repeating the pinned MATLAB controller
generation changed its steady-state and torque tables: the two MATLAB runs
started 0.0164 rpm and 62.8 kN m apart in rotor speed and generator torque,
and their surge traces separated by 6.08 mm at 30 s despite identical
wind, wave, hydro, and mass inputs. The paired Python run now reads each
run's generated controller tables as static case configuration; no saved
MATLAB motion or turbine history enters its solve. The six-coordinate platform,
rotor, and controller solve converges in 14 whole-trajectory passes with
position and velocity iteration residuals below `3e-7` in their respective
units. Against the regenerated source, maximum position differences are
0.477 mm surge, 0.284 mm sway, 0.079 mm heave, `3.41e-5` rad roll,
`2.89e-6` rad pitch, and `1.79e-5` rad yaw. Maximum velocity differences
are 0.0255 mm/s surge, 0.0685 mm/s sway, 0.0443 mm/s heave,
`8.26e-6` rad/s roll, `8.12e-7` rad/s pitch, and `7.29e-6` rad/s yaw.
Rotor speed, azimuth, and generator torque differ by at most `0.000219`
rpm, `0.000299` rad, and 686 N m. The independently coupled three-blade
root loads now have a paired 4%-of-source-peak gate for each of their six
force/moment components; the largest observed component difference is 3.33%
of its source peak. The 30 s gates use 0.6 mm surge,
`1e-5` rad/s yaw speed, and `4e-4` rad azimuth; other trajectory gates
remain as in the 10 s test. The longer
source trace exposed a frame error: its angular velocity is expressed in
world coordinates, while Euler angle rates are different at nonzero pitch
and yaw. The Python platform now maps those rates to world angular velocity
and uses rotated rigid-body inertia and its gyroscopic moment. No wave,
mooring, hydrodynamic, or turbine coefficients were tuned to this source.

A [pinned 60 s MOST developed-sea comparison](https://github.com/cmudrc/wec-sim-python/pull/155)
uses the same 6 m, seed-2 sea and 8 m/s turbulent wind. Its 2,000-frame
TurbSim excerpt supplies 1,249 advected frames through 62.4 s; the first
30 s of the MATLAB trajectory exactly reproduce the separate 30 s source.
Over 6,001 samples, independently synthesized wave elevation and excitation
differ by at most `2.65e-13` m and `3.25e-5` N or N m. The independent
six-coordinate Python solve converges in 19 passes, with final position and
velocity iteration residuals below `1e-6`. Maximum position differences are
0.513 mm surge, 0.632 mm sway, 0.079 mm heave, `5.26e-5` rad roll,
`6.61e-6` rad pitch, and `2.50e-5` rad yaw. Maximum translational velocity
difference is 0.0845 mm/s; maximum angular velocity difference is
`8.26e-6` rad/s. Rotor speed, azimuth, and generator torque differ by at
most `0.000241` rpm, `0.000347` rad, and 686 N m. The largest blade-root
force/moment component error is 3.32% of its source peak, below the 4% gate.
The 60 s position gates are 0.75 mm sway, `6e-5` rad roll, and `8e-6` rad
pitch; all other 30 s trajectory and turbine gates remain in force.

On the 60 s MATLAB motion, Python's restoring, drag, radiation, mooring, and
tower force laws reproduce their logged source channels. Maximum mooring
force and moment errors are `1.54e-5` N and `4.47e-4` N m; the largest tower
moment error is 0.033 N m. Advancing Python's platform with the saved MATLAB
turbine loads still gives 0.630 mm maximum sway difference, nearly the same
as the independently coupled run. A half-step platform diagnostic with
interpolated turbine histories gives 0.674 mm; that interpolation prevents a
strict time-step convergence claim. These comparisons motivate the explicit
60 s gates without changing any hydrodynamic, mooring, or turbine coefficient.
A [pinned published 1,000 s MOST source run](https://github.com/cmudrc/wec-sim-python/actions/runs/38066481378)
uses the full 20,750-frame TurbSim input and the 4 m, seed-1 JONSWAP sea.
MATLAB saves 100,001 wave samples at the 0.01 s solver step and 10,001 body
samples at the published 0.1 s output step. Python's wave elevation differs
by at most `4.84e-13` m, and all six excitation components sampled on the
body grid by at most `2.13e-5` N or N m (gates `1e-12` m and `1e-4` N or N m).
The same-run generated controller reproduces the initial rotor speed and
generator torque to rounding. Its saved restoring, viscous, nonlinear static
mooring, and tower-reaction laws agree across all 10,001 body states within
`5.97e-8`, `9.32e-10`, `4.87e-4`, and `0.0421` N or N m respectively, below
their paired force/moment gates. The BTS SHA-256 is
`060bf3f6bfbebe1dd38f260765611d7d175c6822e9eb3e27e5360ea851df4095`.
A subsequent [fresh full-duration run](https://github.com/cmudrc/wec-sim-python/actions/runs/38071234135)
saved identical wave and excitation histories and the same BTS hash, but its
generated controller had 19 torque-table entries versus 20 in the earlier
artifact. Initial rotor speed and torque changed by `0.01641` rpm and
62.83 kN m; the two MATLAB body trajectories differed by up to 6.44 mm
surge and 3.28 mm sway. Python therefore reads each source run's generated
controller tables as case inputs rather than treating its trajectory as
repeatable across controller regeneration. The full-duration force-law gates
above reconstruct the earlier saved trajectory. A full-history fixed-point
diagnostic over 100 s
was stopped after five passes because its position and velocity residuals
were still `5.24` m and `1.28` m/s. A causal turbine/platform stepper avoids
those repeated full-history passes. The fresh-source [100 s gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38069378327)
passes all six body positions and velocities, rotor speed, azimuth, generator
torque, and three-blade root loads. A local independent 100,001-step solve
against the pinned full 1,000 s source has maximum position differences of
1.74 mm surge, 9.13 mm sway, 0.159 mm heave, `0.000195` rad roll,
`0.0000446` rad pitch, and `0.000359` rad yaw. The largest translational speed
difference is 0.628 mm/s; rotor speed, azimuth, and generator torque differ
by at most `0.00101` rpm, `0.00172` rad, and 4.10 kN m. Each blade-root load
component differs by at most 1.77% of its source peak. The sway error rises
from 0.698 mm over the first 100 s to 9.13 mm over the full 1,000 s, while
the source's peak sway excursion from its initial state is 2.83 m. The
full-duration gate applies explicit per-axis position and velocity limits,
plus rotor, torque, and blade-load limits; the [fresh-source CI run](https://github.com/cmudrc/wec-sim-python/actions/runs/38071234135)
passed all 10,001 saved body samples using that run's generated controller.
The local maxima above refer to the earlier controller realization. This
comparison covers the published 8 m/s turbulent-wind controller and 4 m,
seed-1 sea, not other MOST regimes.

A [pinned 12 m/s MOST constant-wind run](https://github.com/cmudrc/wec-sim-python/pull/142)
changes the published wind-class option and initial steady-state wind speed
but retains the IEA 15 MW turbine, VolturnUS hydro, and 4 m, seed-1 JONSWAP
sea. MATLAB logs `BladePitch` after a Radians to Degrees block in
`MOST_Lib.slx`; Python's controller and BEM use radians. On the saved MATLAB
rotor-speed path, Python reconstructs generator torque within 3.50 N and
converted pitch within `3.2e-7` rad. The independent Python coupling uses
only the run's generated controller tables and the pinned physical inputs.
Over all 1,001 samples and ten coupled passes, maximum platform position and
velocity differences are 0.306 mm and 0.102 mm/s; rotor speed, azimuth,
generator torque, and physical blade pitch differ by at most `0.000375` rpm,
`0.000140` rad, 1.65 kN m, and `1.08e-5` rad. Each blade-root force/moment
component differs by at most 0.084% of its source peak, against a 0.2% gate.
The new causal one-pass solver is independently checked against the same
source and identical per-channel gates. Its local maximum position and speed
differences are 0.305 mm and 0.101 mm/s; rotor speed, azimuth, torque, and
blade pitch differ by at most `0.000389` rpm, `0.000151` rad, 1.71 kN m, and
`7.19e-6` rad. Its largest blade-root component error is 0.085% of the source
peak. The [fresh-source CI gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38075809652)
passed both coupling methods against the same generated controller.
The case starts with a nonzero steady-state pitch setting and exercises pitch
control, but rotor speed subsequently falls and pitch returns to zero. It
covers ten seconds of the 20 s wave ramp, not sustained above-rated operation
or a different published turbulent-wind case.

A [derived 16 m/s, 30 s constant-wind pair](https://github.com/cmudrc/wec-sim-python/pull/164)
keeps the same 4 m, seed-1 JONSWAP sea, platform, and IEA 15 MW turbine while
changing only the wind speed and duration. MATLAB blade pitch remains above
`0.193` rad from 20–30 s, after the wave ramp. On all 3,001 samples, the
independent causal Python solve differs by at most `2.17` mm in platform
position, `0.765` mm/s in platform speed, `0.00255` rpm in rotor speed,
`0.00394` rad in azimuth, `0.601` N m in generator torque, and `0.000607` rad
in blade pitch. The worst blade-root component differs by `0.433%` of its
source peak. Replaying the controller on MATLAB's output-sampled rotor speed
differs in pitch by at most `0.000205` rad, within the `0.0003` rad replay
gate. The longer coupled trajectory has separate explicit per-axis and turbine
gates in `tests/test_most_constant_wind_source.py`; the 10 s gates remain
unchanged. MATLAB's platform travels `17.92` m in surge by 30 s, so this is a
numerical comparison and does not establish physical validity of the source's
small-motion hydrodynamics at that displacement. Turbulent above-rated wind
remains unpaired.

The WaveStar NMPC solver-step audit in
`tools/wavestar_nmpc_step_audit.py` compares the first 30 s of the pinned
published `ode8`/0.05 s source run with independent Python closed loops at
0.001 and 0.0005 s plant steps, both retaining the 0.05 s control step.
Halving the Python step changes pitch by less than `1e-12` rad before NMPC
starts at 15 s and at most `5.82e-5` rad through 30 s. Over the same 15–30 s
window, the 0.001 s Python trajectory differs from the published MATLAB pitch
by up to `0.02877` rad; torque-command differences are `0.00743` N m between
Python steps versus `2.874` N m against MATLAB. The source/Python pitch gap
already reaches `0.00342` rad before any controller acts. Along with the
derived fine-step MATLAB comparison above, this argues against Python plant
step refinement as the remedy for the published coarse-step trajectory gap.
It does not isolate every difference in the active closed loop or establish
full NMPC parity. The physical Python plant and PTO laws remain unchanged.

A [derived fine-step WaveStar NMPC activation run](https://github.com/cmudrc/wec-sim-python/pull/143)
extends the 0.001 s `ode4` source through 15.2 s with the published sea,
geometry, 10 s AR-predictor start, and 15 s NMPC start. The source's
controller sample follows `simu.dt`, so this experiment changes both the
plant and controller steps from the published 0.05 s setting. Direct MATLAB
`ar` and `forecast` calls on the saved 180-sample excitation-estimate window
reproduce the logged 40-step forecasts at six selected active times. The
fit's regression matrix has condition number `6.74e13` at 15 s; the logged
forecast reaches `60.1` kN m within its 40 ms horizon while the estimated
moment stays below `5.35` N m. On saved source inputs, Python reconstructs
the five observer states within `2.8e-12` and the NMPC command within
`4.5e-6` N m. Python's current AR least-squares cutoff discards three
numerically weak directions in that ill-conditioned window; retaining them
improves this derived forecast replay but fails the existing published
composed-controller gate, so the production predictor remains unchanged.
The independent Python run differs over the first 0.2 s of active control by
at most `0.000190` rad pitch, `0.0155` rad/s pitch speed, `5.68` N m torque
request, and `27.8` N PTO force. These are diagnostics, not NMPC parity gates.
The published 0.05 s controller and full closed-loop trajectory remain open.

For the floating OWC, the nine-line MoorDyn replay reconstructs the published
coupling pose and velocity from the floater's center state and its body-local
attachment. It then advances the pinned native MoorDyn library for all 50,001
samples without reading MATLAB's force or line-tension history as input. Each
of six mooring-force components and five named fairlead tensions has a maximum
error below `0.1%` of its source peak and RMS error below `0.01%` of its peak.
The published `lines.txt` repeats the `FairTen3` output label; the MATLAB
struct retains one named channel, so the Python gate aligns channels by name.
The floater motion remains prescribed for this replay. Live feedback to the
seven-coordinate body solve is validated separately below.

The public `solve_floating_owc` now advances the seven-coordinate floater and
water-column state, nine-line native MoorDyn model, chamber pressure, and
turbine speed together from published initial conditions. It uses the pinned
HDF5 regular-wave excitation and body-local mooring attachment, with no saved
MATLAB force, motion, pressure, or speed history used after initialization.
On the local ARM64 500 s run, maximum floater translation error is `0.00285` m,
column translation `0.00065` m, rotation `6.34e-5` rad, pressure `5.48` Pa,
rotor speed `0.237` rad/s, turbine load power `10.93` W, and pneumatic power
`42.93` W over 50,001 samples. Explicit CI gates are `0.01`/`0.005` m for
floater/column position, `0.0002` rad rotation, `15` Pa pressure, `0.5` rad/s
rotor speed, and `25`/`75` W for load/pneumatic power; all six live mooring
force components and five named fairlead tensions must remain within `1%` of
their respective MATLAB peaks. This verifies the published regular-wave,
fixed-frequency-hydrodynamic case, not other seas, radiation-memory settings,
or arbitrary floating OWC designs.

A derived 150 s run keeps that same published sea, geometry, hydrodynamics,
mooring, and air train while setting the slider PTO to `15,000` N/m stiffness
and `60,000` N s/m damping. MATLAB's logged PTO force equals
`-stiffness × stroke - damping × speed` within `3.5e-9` N. The public
`solve_floating_owc` advances the altered system independently, and the same
paired gate now runs through `WEC.floating_owc` and `WEC.run` with a named
floater, column, body-local mooring point, and PTO. The local
maximum differences are `0.00222` m floater position, `0.00086` m column
position, `0.00258` m PTO stroke, `100` N PTO force, `438` W mechanical PTO
power, `8.43` Pa chamber pressure, `0.296` rad/s rotor speed, and `21.5` W
turbine load power over 15,001 samples. The paired CI gates allow `0.007` m
floater position, `0.003` m column position, `0.008` m stroke, `400` N PTO
force, `1,500` W mechanical power, `20` Pa pressure, `0.8` rad/s rotor speed,
and `50` W turbine power; six mooring-force components must stay within
`1.5%` of their respective source peaks. This checks configurable PTO
stiffness/damping in a physically distinct case without fitting source
motion or forces as inputs.

The published RM3 `End_Stops` trajectory has a separate refined-step gate;
the pinned 0.1 s output is not treated as a converged motion reference.

The `MATLAB reference model baselines` workflow runs the two canonical core
examples, all five Sphere free-decay cases, RM3 body-to-body Cases 1–4,
the RM3 MooringMatrix case when selected by workflow dispatch,
both RM3 PTO extension free decays, four RM3 radiation-force baselines (three
paired to Python),
all eight physical conditions from RM3 Multiple Condition Runs Option 1,
three actual `wecSimMCR` runs, each with eight saved body and PTO traces,
the published Sphere passive, reactive PI, declutching, and latching
controllers, a configured Sphere PTO case derived from the passive input,
and the published OSWEC `Nonhydro_Body` case. An earlier
[nine-job run](https://github.com/cmudrc/wec-sim-python/actions/runs/37641710801)
passed against fresh MATLAB R2025b outputs. The Option 1 MCR input was also
expanded into eight scalar simulations. The three driver baselines run
MATLAB's MCR orchestration and postprocessing directly. The pinned Option 2
workbook and Option 3 MAT file select the same ordered cases as Option 1.
Python's MCR interface pairs its case trajectories, PTO signals, average
powers, and power-matrix placement to each driver's actual output.
It records time, position, velocity, total force, and excitation force for
each body as CSV artifacts. New runs also record position, velocity, internal
mechanics force, and power for each PTO. The Applications source is pinned to
[`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/d53d4d4c9eda2581f04204f5d394a6ef84bb099e).
The first run [passed all three model jobs](https://github.com/cmudrc/wec-sim-python/actions/runs/37471757955),
producing nine finite body trajectories with 4,001 samples each (two RM3,
two OSWEC, and five Sphere). The `1m-ME` and `1m` Sphere trajectories were
identical in the saved motion and force signals: the Morison element in that
published case has nonzero coefficients only in x while the free-decay motion
is in heave.
The renewed [model baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/37472556464)
also generated the Sphere HDF5 file with current MATLAB BEMIO and verified
Python preprocessing of it. `wecsim/linearHeave.py` uses that
preprocessing, the heave restoring coefficient, infinite-frequency added
mass, and radiation impulse-response kernel. It integrates the resulting
linear convolution equation with a fixed 0.01 s trapezoidal step. The five
Sphere trajectory comparisons [passed in the same MATLAB job](https://github.com/cmudrc/wec-sim-python/actions/runs/37476318082)
so that the source, HDF5, time grid, and output remain paired.
The [current reference-model run](https://github.com/cmudrc/wec-sim-python/actions/runs/37483219008)
passed all RM3, OSWEC, and Sphere jobs, including Python preprocessing of
the current RM3 and OSWEC HDF5 inputs for both bodies and the focused dynamics
comparisons described above.
The RM3 baseline also checks a reduced Python heave model and a coupled
surge/heave/pitch model. The OSWEC baseline sets the published case's
`waves.phaseSeed` to 1 so repeat runs use one reproducible realization,
saves its realized phase matrix, and checks Python PM binning, wave elevation,
directional excitation, and hinged-pitch motion against that same realization.
The Python wave generator uses NumPy seeds by default. The public
`PMWave(seed=..., phase_generator="matlab")` reproduces pinned Threefry
substreams 1–3 in the paired OSWEC continuous-heading runs and substream 5
in the three-heading Morison run, without a saved phase file. The
case-driven runner shares a generalized dynamics engine across the validated
layouts. A PM phase CSV can replay the actual MATLAB realization.
The RM3 body-to-body job generates the Applications HDF5 input with the pinned
MATLAB BEMIO code, then compares the Python solver with coupling off and on
against paired MATLAB Cases 1 and 2. The CLI selects these hydrodynamic modes
with `rm3` and `rm3 --b2b`; neither command executes an arbitrary application
input file.
The [four-job reference-model run](https://github.com/cmudrc/wec-sim-python/actions/runs/37486493413)
passed RM3, RM3 body-to-body, OSWEC, and Sphere, including the Sphere CLI
smoke test against the MATLAB-generated HDF5 input.
The [case-driven dynamics run](https://github.com/cmudrc/wec-sim-python/actions/runs/37491350883)
passed all four jobs again using JSON cases through the main runner. It also
compared the mapped linear-coordinate RM3 heave and Sphere free-decay cases
against the paired MATLAB trajectories.
Optional linear PTO equilibrium offsets and scalar pretension use the same
generalized dynamics engine. Their nonzero-force behavior has analytical and
case-level tests; the paired MATLAB reference cases use zero offsets.
The configurable small-motion device path also accepts named body motions and
body-local or fixed-world PTO endpoints. Its attachment geometry, projected
stroke, generalized forces, and damping power have analytical and case-level
checks. The configured Sphere case pairs a nonzero spring and extra damper
with MATLAB output while specifying a shifted body-local PTO point. Because
that body only heaves, the shifted point does not affect its trajectory;
attachment-location dynamics remain unpaired.
`wecsim.WEC` provides a Python builder for this same validated path and
returns named NumPy body, coordinate, and PTO histories. The JSON case runner
remains available for saved cases; the Python builder covers the mapped
`linear_subspace`, one-body `floating_gbm`, two-body RM3-style
`floating_joint`, two-body coupled `floating_owc`, and one- or two-body
OSWEC `fixed_hinge` layouts. The
floating joint uses a relative-heave PTO; the fixed hinge accepts a torsional
PTO at the published world-axis locations. Neither accepts arbitrary
attachment geometry.
`python -m wecsim CASE.json --output motion.csv` runs a
supported dynamics configuration. The case declares wave, body, constraint,
PTO, and time settings; the result includes all six body position and
velocity coordinates, applicable wave and PTO signals, and a JSON record of
the case and HDF5 hashes, NumPy version, and Git state. The older
`referenceRunner` remains as a preset CLI for the three model families.
Neither command parses an arbitrary WEC-Sim or Simscape input file.
Expanded application cases remain inventoried. The manually dispatched
`MATLAB reference application regression`
workflow can run the upstream test suites for all 18 application folders
containing the 46 explicit RM3, OSWEC, and Sphere cases. This checks the
pinned MATLAB reference for those cases, but there is no Python time-series
comparison yet.

The [first full application sweep](https://github.com/cmudrc/wec-sim-python/actions/runs/37472556454)
recorded 42 passed methods, two failed Passive Yaw assertions, two MoorDyn
methods filtered by upstream's CI assumption, and an empty Multiple Wave
Spectra suite. Both Variable Hydro methods passed under MATLAB R2025b. A
targeted direct run subsequently passed Multiple Wave Spectra, bringing
the exercised passing checks to 43.

The first application sweep revealed upstream test-suite limitations. In
`Passive_Yaw`, two stored irregular-yaw regression checks fail; the same two
checks fail in [the Applications repository's own CI run](https://github.com/WEC-Sim/WEC-Sim_Applications/actions/runs/36746373598)
on both Linux and Windows under MATLAB R2024b. MoorDyn tests in `Mooring` and
`Paraview_Visualization` explicitly skip themselves on GitHub CI. Forcing the
MoorDyn test initially exposed a native library mismatch (`GLIBCXX_3.4.32`
unavailable from MATLAB's bundled `libstdc++`). Staging the pinned Linux
MoorDyn binaries and preloading the runner's compatible C++ runtime resolved
that reference-environment issue. Targeted R2025b runs now pass the published
[RM3 MoorDyn](https://github.com/cmudrc/wec-sim-python/actions/runs/37884438725)
and [RM3 MoorDyn ParaView](https://github.com/cmudrc/wec-sim-python/actions/runs/37885063313)
tests. The independent Python RM3 body/PTO/mooring trajectory is paired above;
selected published RM3 wave surfaces are paired above; body, mooring, and
full application ParaView histories remain unpaired.

The `Multiple_Wave_Spectra` test class is excluded by MATLAB because its
class name does not match its filename. Our harness generates its OSWEC HDF5
file with BEMIO, runs the input file directly, and requires body output
instead of accepting an empty test suite as success. The
[direct case passed](https://github.com/cmudrc/wec-sim-python/actions/runs/37476616198).
The remaining Passive Yaw findings are MATLAB-reference gaps, not evidence
of Python dynamics agreement or disagreement.
Separate direct runs of the pinned [published irregular passive-yaw input](https://github.com/cmudrc/wec-sim-python/actions/runs/37718464539)
and its [0° heading-update control](https://github.com/cmudrc/wec-sim-python/actions/runs/37720514831)
completed and generated the paired data summarized above. Those runs do not
change the upstream regression assertions.

## Known differences and next reference case

Current MATLAB WEC-Sim changed its PM/JS spectra, seeded phase generator,
object properties, and some wave inputs since the original Python port.
The general `WaveClass` now has direct paired Traditional and EqualEnergy
PM/JS checks using its explicit MATLAB phase generator; the focused OSWEC PM and other
application-specific irregular-wave paths above have separate paired gates.
The historical fixtures alone do not establish parity for remaining wave modes
or all MATLAB random streams and wave modes; the verified Threefry option
applies to the pinned `waveClass` phase path and the listed seeds.
RM3 state-space Cases 5 and 6 remain a radiation-fit investigation: a
passive fit would need its own evidence against source BEM data and
convolution trajectories before it could count as physically validated.
The [pinned Applications `End_Stops` input](https://github.com/WEC-Sim/WEC-Sim_Applications/blob/d53d4d4c9eda2581f04204f5d394a6ef84bb099e/End_Stops/wecSimInputFile.m) sets
`upperLimitTransitionRegion` and `lowerLimitTransitionRegion` to 0.5 m, but
the [pinned `ptoClass.hardStops` defaults](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/ptoClass.m) and the referenced translational PTO
Simulink block read the fields ending in `TransitionRegionWidth`. The input
therefore leaves the effective widths at their 1e-4 m defaults. The
[fresh source export](https://github.com/cmudrc/wec-sim-python/actions/runs/37689167739)
confirms these effective settings and supplies body/PTO trajectories,
hydrodynamic force components, and original/applied mass matrices. The PTO
stroke reaches −0.684 to +0.687 m versus about ±0.865 m in the published
no-stop case, so contact materially changes the trajectory. A refined 400 s
source run records two samples inside the effective transition width; the
Python `LinearHardStops` smoothing law reconstructs its logged PTO force to
within `2.4e-7` N. The upstream
[test](https://github.com/WEC-Sim/WEC-Sim_Applications/blob/d53d4d4c9eda2581f04204f5d394a6ef84bb099e/End_Stops/TestEndStops.m)
only checks that `wecSim` runs.

The source regular-wave case applies split added mass with delayed
acceleration feedback. Reconstructing its translational added-mass force from
the two preceding 0.1 s output accelerations agrees to numerical precision
until 84.4 s; 78 later output samples differ by more than 0.001 N, with a
maximum residual of 17.9 MN. Those intervals are consistent with internal
variable solver steps near impacts, whose acceleration history is absent from
the exported 0.1 s samples. The [0.05 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37691660845)
and [0.025 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37692499385)
and [0.0125 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37693255891)
have no such saved-step residual above 0.001 N through 120 s. PTO stroke
changes by up to 20.7 mm from 0.05 to 0.025 s and 2.79 mm from 0.025 to
0.0125 s, versus 149.4 mm from 0.1 to 0.05 s. Ordinary PTO damper energy
over 120 s falls from 11.176 MJ in the published 0.1 s run to 10.518 MJ at
0.0125 s.

Over the full 400 s, the two refined MATLAB runs differ by at most 2.79 mm
PTO stroke, 18.8 mm/s PTO speed, and 242.4 kN PTO force. Their ordinary
damper energies are 55.062 and 54.823 MJ, a 0.437% difference. Neither
refined run has a saved-step added-mass recurrence residual over 0.001 N;
the published 0.1 s run has 78 such samples and absorbs 60.941 MJ. The
source convergence gate also checks both bodies' active positions and
velocities. These are numerical step checks, not altered WEC settings.

The Python hard-stop model instead applies added mass implicitly and resolves
contact forces between output samples. Its 0.025 and 0.0125 s output settings
change PTO stroke by less than 1 µm over 400 s. This physical dynamics
path is available through the RM3 solver and case runner; it does not embed
the source's delayed acceleration loop. Against the 0.0125 s MATLAB run, the
Python body's largest surge and heave position errors are 1.08 and 0.936 mm;
the largest pitch error is `6.74e-5` rad. Ordinary damper energy is 54.747 MJ,
0.0762 MJ (0.139%) below the refined source. General nonlinear stop damping
and arbitrary WEC layouts remain outside this RM3 validation.
Further dynamics targets include other published OSWEC and Sphere
configurations. The Sphere free-decay cases already have direct Python motion
comparisons. A comparison must
record both code revisions, the HDF5 input, time step, outputs, and numerical
tolerances.

The inherited `paraviewClass.py` was invalid Python mixed with unfinished
MATLAB code. It now calculates wave surfaces and writes source-paired wave,
triangular-body, and mooring-line VTP files. The actual published RM3 wave
meshes are paired at three times, and the OSWEC nonlinear flap pressures are
paired separately. Other application-specific geometry and pressure histories,
actual line histories, and full visualization collections remain open.
The older object tests import duplicate copies of classes inside test folders.
The parity tests here import the production files instead.
Six inherited files under `tests/test_simulink` still contain unfinished
MATLAB-like Python and do not parse; they are outside the collected pytest
suite. The `wecsim` package compiles, and the scoped production parity suite runs.
