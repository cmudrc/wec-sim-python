# wec-sim-python

> **cmudrc fork status:** This is an active parity effort, not yet a complete
> wave energy converter simulator. Selected Sphere, RM3, OSWEC, and barge
> cases have paired MATLAB trajectory checks; general WEC-Sim dynamics and
> several published applications remain unsupported.
> See [PARITY.md](PARITY.md) for verified behavior, current
> MATLAB reference revision, and the remaining work.

Install from a clone with Python 3.12:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -c 'from wecsim import WEC; print(WEC)'
```

For a regular install without cloning, use
`python -m pip install git+https://github.com/cmudrc/wec-sim-python.git`.
Both methods install the `wecsim` import and the `wecsim` and
`wecsim-reference` commands.

For development, run the production-code parity checks with:

```sh
.venv/bin/python -m pip install pytest
.venv/bin/python -m compileall -q wecsim
.venv/bin/python -m pytest -q tests/test_wave_parity.py tests/test_body_io.py tests/test_oswec_standalone.py tests/test_rm3_standalone.py tests/test_reference_cli.py tests/test_general_dynamics.py tests/test_case_dynamics.py tests/test_pto_connections.py tests/test_python_api.py
```

**WEC-Sim-Python** is Sungjun Won's Python port of
[WEC-Sim](https://github.com/WEC-Sim/WEC-Sim), the MATLAB/Simulink wave energy
converter simulator. This fork is developing and checking the Python code
against MATLAB WEC-Sim while preserving the original author's work.

The installable code is in `wecsim/`, and runnable cases and bundled RM3
inputs are in `examples/`. Sungjun Won's early sketches remain available in
Git history.

## Goal of WEC-Sim-Python
**WEC-Sim-Python** aims to help researchers, start-up companies, and enthusiasts without access to MATLAB in order to use the open-source code provided by NREL and Sandia lab. Also, with growing research in the field of machine learning, **WEC-Sim-Python** could be more convenient for those who develop machine learning projects utilizing Python.

## Current status

Wave generation, RM3 and OSWEC hydrodynamic input, force preprocessing, and
the supported device dynamics have paired MATLAB checks. The published
generalized-body-mode barge now has a coupled rigid and flexible regular-wave
runner for its floating three-DOF joint. Other GBM layouts and wave options
remain unsupported; see [PARITY.md](PARITY.md). The Python API is
the primary way to configure a supported device. It constructs bodies,
named motions, attachment points, PTOs, and waves as Python objects, then
returns NumPy arrays directly. Import it as `wecsim` after installation:

```python
from wecsim import NoWave, WEC, WorldPoint

wec = WEC("Heaving float")
float_body = wec.body("float", "path/to/hydro.h5")
wec.coordinate("heave", float_body.move("heave"))
wec.pto("main", float_body.at(0, 0, 0), WorldPoint(0, 0, 10),
        damping=1_200_000)
result = wec.run(NoWave(), dt=0.1, end_time=20,
                 radiation_memory=15,
                 initial_coordinate={"heave": 1})
heave = result.bodies["float"].position[:, 2]
pto_force = result.ptos["main"].force
```

The [full Python example](examples/configurable_rm3_pto.py) defines a
two-body device with a shared pitch pivot and body-local PTO points. Run it
with `python -m examples.configurable_rm3_pto`. Relative HDF5 paths in
the Python API resolve from the current directory unless `base_dir` is passed
to `wec.run`. Body order must match the HDF5 hydrodynamic body order. The
Python builder covers `linear_subspace`, the single-body `floating_gbm`
regular-wave layout, the paired two-body `floating_joint`, and the OSWEC
`fixed_hinge` layout. Mapped
coordinates use small-motion kinematics and fixed-axis PTOs. The floating
joint uses its own pitched-slider geometry and a relative-heave PTO; arbitrary
PTO attachment points are not part of that reduced layout.

The MOST application uses a TurbSim full-field wind input. Its checked-in
`.bts` file can be read directly with `wecsim.read_turbsim_bts(path)`; the
result has `velocity[time, component, y, z]`, `y`, `z`, `dt`, and `hub_height`.
`wecsim.MostWindField(wind).at_index(i)` evaluates the ten published X planes
at one output time as `[component, x, y, z]`, without storing the full
five-dimensional shifted wind array. This covers the MOST wind input and its
frozen-turbulence time shift. `MostWindField(wind).sampler(i)` returns the
trilinear world-position wind callable for that output time, suitable for
`MostBEM.loads(hub_state, pitch, sampler)`.
`MostWindField(wind).sampler_at(time)` also interpolates between advection
frames at the simulation time in seconds.
For MOST's constant-wind option, use `wecsim.MostConstantWind(12,
end_time=10)` in place of `MostWindField`. Its default direction is world
positive X; `direction` and `domain_sizes` can be configured to match a
constant-wind MATLAB input.
`wecsim.MostStaticMooring().force(pose)` evaluates
the published three-line static catenary load and each fairlead's horizontal
and vertical tension for a six-coordinate platform pose. For the published
IEA 15 MW turbine, `wecsim.MostBEM.from_iea15mw(blade_data_directory)` reads
its raw AeroDyn blade and airfoil tables; `model.loads(hub_state, pitch, wind)`
returns six root loads per blade for a 14-state hub, blade pitch in radians,
and either a wind three-vector or a position-to-wind callable.
`wecsim.MostBaselineController.iea15mw().simulate(time, rotor_speed)` returns
the published baseline generator torque and blade pitch for a prescribed
rotor-speed history. `wecsim.MostRotor.from_iea15mw(blade_data_directory)`
advances the published rotor speed, azimuth, generator torque, and blade pitch
from a supplied six-DOF platform position and velocity history and a
`MostWindField` or `MostConstantWind`:

```python
from wecsim import MostRotor

rotor = MostRotor.from_iea15mw(blade_data_directory)
result = rotor.simulate(time, platform_position, platform_velocity, wind_field)
```

To reproduce a particular MATLAB MOST run, pass its generated control tables
as case configuration. The source steady-state optimization can produce
slightly different tables on separate runs:

```python
from wecsim import MostBaselineController

controller = MostBaselineController.from_matlab_files(
    control_mat, steady_states_mat, wind_speed=8.0,
)
rotor = MostRotor.from_iea15mw(blade_data_directory, controller=controller)
```

Rotor speed in `result` is rad/s. The MATLAB MOST output reports rpm and logs
blade pitch in degrees after a conversion block; Python blade pitch is in
radians. This runner takes platform motion as input; turbine mass and
aerodynamic loads do not yet feed back into independently advanced platform
motion.
Its `blade_root_load` array contains the BEM loads at the advanced rotor
states, in the preconed blade frames.

`MostTowerReaction.from_iea15mw(properties_file).evaluate(position, velocity,
acceleration, rotor_speed, azimuth, generator_torque, blade_root_load)`
computes the six-component force and moment transferred from this turbine to
the platform. It uses the published tower, nacelle, hub, and blade mass
properties, gravity, rotor inertia, and blade-root aerodynamic loads.
`blade_root_load` is an N×6×3 array in the preconed blade frames; the returned
N×6 wrench is in the rotating platform frame about the tower base. This is a
force calculation on supplied platform and rotor states, not a coupled solve.

`MostPlatformHydrodynamics.from_volturnus(h5_file, mass_properties_file)`
reads the published VolturnUS BEMIO coefficients and platform mass properties.
Its `restoring_force(position)`, `drag_force(velocity)`, and
`radiation_force(velocity, dt)` evaluate the six-component forces on a
platform trajectory, using WEC-Sim's logged resisting-force signs. Position
and velocity are arrays with six columns in surge/sway/heave/roll/pitch/yaw
order. `platform.simulate(time, wave_excitation, tower_base_load)` advances the
published dominant surge, heave, and pitch coordinates from these Python
forces and the live nonlinear-static mooring. `wave_excitation` is the
six-component output of Python's JONSWAP/HDF5 synthesis;
`tower_base_load` is a six-component history in the platform frame. This
one-way solve takes tower reaction as input and does not yet feed Python
turbine loads back to the platform.

For a coupled trajectory, use `MostCoupled(platform, rotor, tower)`
with a Python-generated six-component wave-excitation history and the MOST
wind field:

```python
from wecsim import MostCoupled

coupled = MostCoupled(platform, rotor, tower)
result = coupled.simulate(time, wave.excitation_force, wind_field)
platform_motion = result.platform
rotor_motion = result.rotor
```

The runner advances its own rotor, BEM loads, controller, tower reaction, and
all six platform coordinates until the two trajectories agree within the
reported `position_residual` and `velocity_residual`. Pass
`full_six_dof=False` for the earlier surge/heave/pitch reduction. The paired
six-coordinate gates cover two pinned 10 s MOST seas, 30 and 60 s developed
seas, a derived 10 s constant-wind case at 12 m/s, and a 30 s constant-wind
case at 16 m/s that retains active blade pitch after the wave ramp.

For longer histories, `coupled.simulate_causal(time, wave.excitation_force,
wind_field)` advances the turbine and platform once per 0.01 s step using a
predicted platform pose. Its reported residuals are the maximum correction
between that prediction and the solved platform state. The 100 s prefix of
the published MOST case passes a fresh MATLAB source gate. An independent
Python run over the full 1,000 s differs by at most 9.13 mm in body position,
0.00101 rpm in rotor speed, and 1.77% of each blade-load component's source
peak. The [fresh 1,000 s source gate](https://github.com/cmudrc/wec-sim-python/actions/runs/38071234135)
passes all saved body samples with its own generated controller tables; see
`PARITY.md` for the per-axis limits, controller-run variability, and remaining
case coverage.

For the published RM3 floating joint, configure the two bodies in HDF5 order:

```python
from wecsim import RegularCICWave, WEC

wec = WEC("RM3 floating joint")
float_body = wec.body("float", "path/to/rm3.h5",
                      inertia=(0, 21_306_090.66, 0))
spar = wec.body("spar", "path/to/rm3.h5",
                inertia=(0, 94_407_091.24, 0))
wec.floating_joint(float_body, spar, damping=1_200_000)
result = wec.run(RegularCICWave(2.5, 8), dt=0.1, end_time=400,
                 ramp_time=100, radiation_memory=60)
stroke = result.ptos["relative_heave"].stroke
```

The default uses implicit added mass. For comparisons with the pinned MATLAB
Simulink block, `added_mass_scheme="simulink_delay"` selects its numerical
feedback setting explicitly. `floating_joint` also accepts a joint surge
spring, PTO stiffness and equilibrium, hard stops, and convolution/FIR
radiation where the underlying paired solver supports them. The returned
`absorbed_power` counts the linear damper; spring energy exchange can be
computed from `-force * velocity`. The returned stroke is float heave minus
spar heave, so an initial body offset appears in the first stroke sample.

For the published generalized-body-mode barge, generate its HDF5 with the
Applications `Generalized_Body_Modes/hydroData/bemio.m`, then run:

```python
from wecsim import RegularWave, WEC

wec = WEC("GBM barge")
barge = wec.body("barge", "path/to/barge.h5",
                 inertia=(6.667e7, 2.167e9, 2.167e9))
wec.floating_gbm(barge)
result = wec.run(RegularWave(height=2, period=8),
                 dt=0.05, end_time=400, ramp_time=100)
surge = result.bodies["barge"].position[:, 0]
mode_displacement = result.flexible_modes["barge"].position
```

The published OWC orifice can be connected to its one flexible body mode:

```python
from math import pi
from wecsim import OrificePTO, PMWave, WEC

wec = WEC("OWC")
body = wec.body("OWC", "path/to/test17a_clean.h5",
                inertia=(99.28, 11.04, 99.2))
orifice = OrificePTO(piston_area=pi * 0.25**2,
                     orifice_area=pi * 0.01**2,
                     discharge_coefficient=0.62, air_density=1.2)
wec.floating_gbm(body, orifice=orifice,
                 heave_linear_damping=100, mode_linear_damping=100,
                 heave_drag_cd=1.2, heave_drag_area=8,
                 pitch_drag_cd=1.2, pitch_drag_area=8)
result = wec.run(PMWave(height=1, period=4, seed=128),
                 dt=0.005, end_time=6, ramp_time=10,
                 radiation_memory=15)
piston = result.ptos["orifice"]
power_watts = piston.absorbed_power
mach_flag = dict(result.raw.extra_outputs)["orifice_compressibility_flag"]
```

The rigid joint has surge, heave, and pitch; the orifice force acts on the
flexible mode and reacts on rigid heave. The settings above reproduce the
published OWC input file. `seed=128` makes a reproducible **Python** sea;
`phase_generator="matlab"` selects the pinned MATLAB Threefry substream.
`compressibility_flag` marks speeds beyond the chosen Mach threshold without
changing the incompressible law. Paired motion is gated only through 6 s;
the Python and MATLAB heave and flexible trajectories diverge later. In the
pinned OWC Simulink model, the piston reaction enters the reported flexible
acceleration but does not enter the flexible state-space integrator. Python
applies that reaction to the flexible motion by default. For a paired source
comparison, pass `orifice_force_path="published_owc"` to `floating_gbm`; this
reproduces the published force routing and keeps the source's omitted piston
force explicit in the case provenance. The published MATLAB air model also
first exceeds its Mach threshold at 10.13 s, so its later trajectory is not a
physical-accuracy target. In the source-compatible setting, reported PTO power
is the source algebraic power signal, not energy removed from the simulated
flexible state.

The separate published `OWC/FloatingOWC` device has a compressible air
chamber. Its pressure state and axial reaction can be configured in Python:

```python
from math import pi
import numpy as np
from wecsim import FloatingOwcChamber, FloatingOwcColumnJoint, FloatingOwcTurbine

joint = FloatingOwcColumnJoint(center_separation=29.445)
floater_pose = np.array([0, 0, -31.945, 0, 0, 0])
floater_velocity = np.zeros(6)
column_pose = joint.column_pose(floater_pose, stroke=0.2)
column_velocity = joint.column_velocity(
    floater_pose, floater_velocity, stroke=0.2, stroke_speed=0.1,
)

area = pi * 5.89**2 / 4
chamber = FloatingOwcChamber(
    area=area, initial_volume=area * 4.5, gamma=1.4,
    ambient_pressure=101325, ambient_density=1.25,
    turbine_diameter=0.75, turbine_kappa=0.775,
)
pressure = 1000.0                 # Pa gauge
water_column_displacement = column_pose[2] - (-2.5)  # m world heave
water_column_speed = column_velocity[2]             # m/s world heave speed
turbine_speed = 150.0            # rad/s
pressure_rate = chamber.pressure_derivative(
    pressure, water_column_displacement, water_column_speed, turbine_speed,
)
column_force = chamber.force_on_column(pressure)
turbine = FloatingOwcTurbine(
    diameter=0.75, inertia=3.06,
    control_coefficient=2e-4, control_exponent=3,
    max_control_torque=216.5,
)
turbine_state = turbine.evaluate(pressure, turbine_speed)
rotor_acceleration = turbine_state.speed_derivative
load_power = turbine_state.load_power
```

Pressure is gauge Pa. The published chamber input is the column's world heave
relative to its equilibrium center, positive when the column rises into the
chamber. The PTO stroke is relative to the floater and can differ substantially
from that world heave. The turbine performance curves are the
published Wells fits; rotor inertia,
control gain, exponent, and torque limit are configurable. ``load_power`` is
the published model's `P_turb = control_torque * turbine_speed` in W. The
chamber and turbine states are paired with the MATLAB application when driven
by its water-column motion. The complete published regular-wave device can
also be advanced with live native MoorDyn, body, chamber, and turbine feedback:

```python
from wecsim import MoorDyn, solve_floating_owc

mooring = MoorDyn("path/to/libmoordyn.so", "path/to/writable/lines.txt")
response = solve_floating_owc(
    "path/to/floatingOWC.h5", mooring,
    wave_height=4.5, wave_period=11.2, ramp_time=50,
    moordyn_point=(0, 0, 31.945),  # body-local point relative to floater CG
    chamber=chamber, turbine=turbine,
    dt=0.01, end_time=500,
)
print(response.column_pose[-1], response.turbine_power[-1])
```

Set `pto_stiffness` and `pto_damping` (N/m and N s/m) in
`solve_floating_owc` to configure the axial PTO. Positive values oppose the
water-column stroke and speed; the response reports axial force, mechanical
power entering the PTO, and damper-only dissipation. Both default to zero, as
in the published application.

The HDF5 and native MoorDyn library are external inputs. The solver defaults
to the published body masses and inertias; wave, body, attachment, chamber,
and turbine parameters can be changed through Python arguments. This solver
uses the published regular-wave fixed-frequency hydrodynamic formulation,
with body-to-body radiation coupling disabled as in the source case. Its
paired trajectory and remaining scope are in [PARITY.md](PARITY.md).

The same coupled dynamics can be configured with the public `WEC` builder:

```python
from wecsim import MoorDyn, RegularWave, WEC

wec = WEC("Floating OWC")
floater = wec.body("floater", "floatingOWC.h5", hydro_body=1,
                   inertia=(1.531e9, 1.531e9, 0.1118e9))
column = wec.body("column", "floatingOWC.h5", hydro_body=2,
                  mass=4_493_450)
wec.floating_owc(
    floater, column,
    moordyn=MoorDyn("path/to/libmoordyn.so", "path/to/writable/lines.txt"),
    moordyn_point=floater.at(0, 0, 31.945),
    column_height=50.69, column_diameter=5.89,
    pto_name="slider", pto_stiffness=15_000, pto_damping=60_000,
)
result = wec.run(RegularWave(4.5, 11.2), dt=0.01,
                 end_time=150, ramp_time=50)
print(result.bodies["column"].position[-1], result.ptos["slider"].force[-1])
```

`result.ptos["slider"].absorbed_power` is damper dissipation. The complete
spring-plus-damper mechanical power, chamber pressure, turbine outputs, and
MoorDyn coupling histories are in `dict(result.raw.extra_outputs)`. The
column inertia defaults to that of a uniform cylinder when omitted.

For a regular-wave declutching PTO, pass
`control=DeclutchingControl(gain=232_020, declutch_time=0.8)` to `wec.pto`
instead of constant damping. Import `DeclutchingControl` from `wecsim`.
The controller disengages at a velocity sign change and reengages after the
configured interval; `minimum_on_time` defaults to 0.2 s. The published
Sphere case is paired against MATLAB at a 0.01 s step. Other geometries and
wave settings are not yet paired.
For the published Sphere latching case, pass
`control=LatchingControl(gain=49_181, latch_damping=37_308_296, latch_time=2.4)`.
The controller applies the larger damping for the timed interval after a
velocity reversal, then returns to its normal gain. It is a finite damping
force, not a rigid lock. `minimum_normal_time` defaults to 0.2 s. The control
settings serialize through `WEC.to_case` for reproducible saved cases.

The JSON runner remains available for saved and reproducible cases. It
accepts simulation, wave, body, constraint, and PTO settings. The included
RM3 example runs with the historical bundled HDF5:

```sh
python -m wecsim examples/rm3.json --output results/rm3.csv
```

The runner saves body position and velocity, wave elevation where applicable,
PTO force or torque, and a JSON provenance record beside the CSV. HDF5 paths
in the case file resolve relative to that file. For a comparison with current
MATLAB WEC-Sim, point the case at the HDF5 file from the pinned reference run.
The provenance record hashes the case, hydrodynamic files, and any replayed
phase CSV.

Supported combinations are:

| Constraint `kind` | Wave `type` | Bodies and PTO | Integration |
| --- | --- | --- | --- |
| `heave` | `none` | One equilibrium-mass body, initial heave displacement, no PTO | Radiation convolution |
| `fixed_hinge` | `pm` or `pm_multi` | Hydrodynamic flap, optional fixed hydrodynamic base, pitch PTO | Directional PM excitation and radiation convolution; `pm_multi` sums independently phased seas |
| `fixed_hinge` | `spectrumImportFullDir` | Hydrodynamic flap, optional fixed hydrodynamic base, pitch PTO | Imported frequency-dependent directional spectrum and radiation convolution |
| `fixed_hinge` | `regular` | One hydrodynamic flap, optional fixed nonhydrodynamic base, pitch PTO | Regular-wave excitation and constant-frequency radiation |
| `fixed_morison` | `pm` | Stationary bodies without HDF5; body-local Cartesian Morison elements; no PTO | Directional irregular-wave velocity, acceleration, and six-component Morison force; optional single-heading current |
| `linear_subspace` | `none` or zero-heading `regular` | One hydrodynamic body with axial body-local Morison elements; pure heave, or surge/heave/pitch in regular waves; no PTO | Radiation convolution or constant-frequency radiation, relative-fluid drag, fluid inertia, and Morison added mass in the acceleration solve; optional regular-wave current for surge/heave/pitch |
| `floating_joint` | `regular` or `regularCIC` | Two equilibrium-mass bodies, pitch inertias, relative-heave PTO; optional `body_to_body` | Constant-frequency radiation, impulse-response convolution, or sampled FIR radiation |
| `floating_joint` | `elevationImport` | Two equilibrium-mass bodies, relative-heave PTO, optional joint surge spring | Imported MAT elevation and radiation convolution |
| `floating_joint` | `none` | Two equilibrium-mass bodies, named initial coordinates and speeds, relative-heave PTO | Radiation convolution for paired free decay; sampled FIR is also available |
| `linear_subspace` | `regular`, `regularCIC`, `pm`, `jonswap`, `spectrumImport`, `elevationImport`, or `none` | Any number of six-DOF hydrodynamic bodies; named motions or 6-by-N maps; optional linear or rotational PTOs; selected mean-drift coefficients for regular waves | Constant-frequency radiation or convolution |

The canonical OSWEC pitch hinge also has a direct Python configuration:

```python
from wecsim import PMWave, WEC, WorldPoint

wec = WEC("OSWEC pitch flap")
flap = wec.body("flap", "oswec.h5", mass=127_000,
                inertia=(1.85e6,) * 3)
base = wec.body("base", "oswec.h5", mass=999, inertia=(999,) * 3)
wec.fixed_hinge(flap, base, location=WorldPoint(0, 0, -10),
                pto_location=WorldPoint(0, 0, -8.9), damping=12_000)
sea = PMWave(2.5, 8, directions=(0, 30, 90),
             spreading=(0.1, 0.2, 0.7), seed=1)
result = wec.run(sea, dt=0.1, end_time=400, ramp_time=100,
                 radiation_memory=30)
pitch = result.coordinates["pitch"].position
torque = result.ptos["hinge"].force
```

The hinge can also take a `fixed_body(...)` base. Both hinge and PTO points
currently lie on the world z axis. PTO angle and speed are in radians and
radians per second; `absorbed_power` reports positive damping dissipation.

The published OSWEC passive-yaw cases use one yaw coordinate and a torsional
PTO. Set `passive_yaw=True` on the moving body to interpolate excitation at
the wave heading relative to its current yaw angle:

```python
from wecsim import RegularWave, WEC, WorldPoint

wec = WEC("Passive-yaw OSWEC")
flap = wec.body("flap", "oswec.h5", mass=12_700,
                inertia=(1.85e6,) * 3, passive_yaw=True)
wec.body("base", "oswec.h5", mass=999, inertia=(999,) * 3)
yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
wec.rotational_pto("hinge", yaw, damping=120_000)
result = wec.run(RegularWave(2.5, 8, direction=10),
                 dt=0.01, end_time=600, ramp_time=100)
angle = result.bodies["flap"].position[:, 5]
absorbed_power = result.ptos["hinge"].absorbed_power
```

The rotational PTO reports angle in radians as `stroke`, angular speed in
rad/s as `velocity`, and torque in N m as `force`. The same device accepts a
Pierson–Moskowitz sea with radiation memory:

```python
from wecsim import PMWave

result = wec.run(PMWave(2.5, 8, direction=10, seed=1,
                        phase_generator="matlab"),
                 dt=0.01, end_time=250, ramp_time=100,
                 radiation_memory=40)
```

`PMWave(..., phase_file="phases.csv")` replays a saved frequency-by-direction
phase matrix. A seed uses NumPy by default; `phase_generator="matlab"` uses
the pinned WEC-Sim Threefry stream, so no phase file is needed for those
seeded cases. Passive yaw supports one pure-yaw
hydrodynamic body, stationary additional bodies, and independent radiation
with full-circle BEM headings. Python interpolates the wave heading
continuously. The published MATLAB regular and irregular passive-yaw inputs
hold heading coefficients until yaw changes by 0.01° and 1°, respectively.
For the published irregular PM case, set `passive_yaw_threshold=1` in the
moving body's `wec.body(...)` call. This reproduces the source's one-degree
coefficient hold and snap to a nearby tabulated BEM heading. A zero threshold
keeps continuous interpolation. The sampled setting is available for PM waves
with one incident direction; small trajectory differences can change its
update sample and accumulate over long runs.
A three-seed paired diagnostic that prescribes only MATLAB's heading-update
schedule recovers tight yaw and PTO-work agreement while Python advances its
own motion. Native long-trajectory parity for the one-degree hold remains open.

For the published `Variable_Hydro/Passive_Yaw` case, set
`passive_yaw=True` and `yaw_heading_bank=np.arange(-30, 30.25, 0.25)` on the
flap. This selects the nearest BEM heading from wave direction relative to
yaw, matching the input's 241 selected headings. For example, the irregular
PM configuration uses the same body, yaw coordinate, and rotational PTO:

```python
import numpy as np
from wecsim import PMWave, WEC, WorldPoint

wec = WEC("variable-yaw OSWEC")
flap = wec.body("flap", "oswec.h5", mass=12_700,
                inertia=(1.85e6,) * 3, passive_yaw=True,
                yaw_heading_bank=np.arange(-30, 30.25, 0.25))
wec.body("base", "oswec.h5", mass=999, inertia=(999,) * 3)
yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
wec.rotational_pto("hinge", yaw, damping=120_000)
result = wec.run(PMWave(2.5, 8, direction=10, seed=1,
                        phase_generator="matlab"),
                 dt=0.01, end_time=600, ramp_time=100,
                 radiation_memory=40)
```

The full 600 s regular-wave yaw, six excitation channels, and PTO work pair
with MATLAB. The same bank
works with a single-direction `PMWave`; a derived 120 s irregular case pairs
the source excitation within `1e-6` N or N m and the independent yaw within
`0.002` rad. A longer MATLAB job completed the published 600 s irregular
case: source-heading selection and excitation pair, while independent Python
yaw differs by up to `0.1103` rad and PTO work by `9.05%`. Prescribing only
MATLAB's heading choices cuts the work difference to `0.0081%`; native
full-length trajectory parity remains open. The source generates a
finer 0.05° set of HDF5 files, but its input reads only the selected 0.25°
headings. The paired baseline generates those same selected files with the
pinned interpolation rule; see [`PARITY.md`](PARITY.md).

The published `Morison_Element/morisonElement` application has a fixed
monopile and tower without HDF5 hydrodynamic bodies. Its nonlinear drag and
fluid-inertia force can be configured with a body-local point:

```python
import numpy as np
from wecsim import PMWave, WEC

wec = WEC("fixed monopile")
monopile = wec.fixed_body("monopile", center_gravity=(0, 0, -15),
                          volume=np.pi * 10**2 * 30,
                          inertia=(1.25e9, 1.25e9, .15e9))
wec.fixed_body("tower", center_gravity=(0, 0, 25),
               mass=1_031_930, inertia=(9.66e8, 9.66e8, .132e8))
wec.morison_element(monopile, point=monopile.at(0, 0, 10),
                    drag_coefficient=(1, 1, 1),
                    added_mass_coefficient=(1, 1, 1),
                    area=(300, 300, np.pi * 10**2 / 4),
                    volume=np.pi * 10**2 * 30,
                    phase_mode="matlab_shared")
sea = PMWave(2, 5, seed=5, directions=(0, 30, 90),
             spreading=(.1, .2, .7), frequency_range=(.001, 10),
             water_depth=30)
result = wec.run(sea, dt=.01, end_time=400, ramp_time=100, rho=1025)
force_and_moment = result.body_forces["monopile"]  # time × 6, N and N m
```

The separate published `Morison_Element/monopile` application uses an HDF5
hydrodynamic monopile and a fixed tower. For that stationary layout, the
public case runner returns wave excitation, both bodies' total forces, and
both fixed-joint reactions:

```python
from wecsim import run_fixed_hydro_monopile
from wecsim.irregularWave import pm_equal_energy_components

hydro = "applications/Morison_Element/hydroData/monopile.h5"
sea = pm_equal_energy_components(
    hydro, significant_height=2, peak_period=5,
    directions=(0,), spreading=(1,), count=500, seed=5,
)
response = run_fixed_hydro_monopile(hydro, sea)
surge_force = response.excitation_force[:, 0]  # N
joint_reactions = response.joint_force  # time × two joints × six components
```

The HDF5 file is generated by the MATLAB application's `hydroData/bemio.m`;
it is not bundled with the Python package. A Python `seed` creates a
reproducible sea. Use `PMWave(seed=5, phase_generator="matlab")` to reproduce
the pinned source sea without a saved phase array. This runner covers the published stationary geometry and
fixed-joint force conventions; it does not solve moving hydrodynamic bodies.

For a fixed Morison body in a **single-heading** PM sea, set a horizontal
current on the wave. Its direction is in world coordinates; `depth` is the
distance below the surface where the profile reaches zero. Current and wave
velocity use the same startup ramp, and current contributes drag but no fluid
acceleration:

```python
from wecsim import Current, PMWave

sea = PMWave(2, 5, directions=(30,), spreading=(1,),
             frequency_range=(.001, 10), water_depth=30,
             current=Current(speed=.8, direction=45,
                             profile="power", depth=30))
result = wec.run(sea, dt=.01, end_time=20, ramp_time=10)
```

`profile` may be `"uniform"`, `"power"` (the one-seventh law), or `"linear"`.
The latter two require `depth`. Active current with multiple wave headings
raises an error because the pinned MATLAB force function adds the current
once per heading.

The pinned MATLAB function uses the first phase column for Morison force at
every heading, while wave elevation uses each heading's own phase.
`phase_mode="matlab_shared"` selects that source behavior; the default
`"directional"` uses each heading's own phase for force. Both modes retain
WEC-Sim's per-heading drag calculation. Supply a saved
three-column phase CSV to reproduce a MATLAB realization exactly.
An axial element on a moving hydrodynamic body can also contribute drag,
fluid inertia, and added mass. For the published Sphere's surge, heave, and
pitch joint in zero-heading regular waves:

```python
from wecsim import Current, RegularWave, WEC

wec = WEC("Sphere with axial Morison element")
sphere = wec.body("sphere", "sphere.h5",
                  inertia=(20907301, 21306090.66, 37085481.11))
for axis in ("surge", "heave", "pitch"):
    wec.coordinate(axis, sphere.move(axis))
wec.morison_element(sphere, point=sphere.at(0, 0, -2),
                    drag_coefficient=(0, 0, 1),
                    added_mass_coefficient=(0, 0, 1),
                    area=(0, 0, 100), volume=20)
result = wec.run(RegularWave(1, 8), dt=.01, end_time=40, ramp_time=10,
                 initial_coordinate={"heave": 1})
morison_wrench = result.body_forces["sphere"]  # time × 6, N and N m
```

The same surge/heave/pitch layout accepts a regular-wave current, for example
`RegularWave(1, 8, current=Current(.8, 45, "power", 30))`. Current loads the
moving element through relative velocity and drag; it adds no fluid
acceleration. The direct source force law is paired for three current profiles.
A derived 40 s Sphere run with current shows a convention gap at nonzero
pitch: the pinned MATLAB Cartesian option-1 function adds no heave force
from horizontal current, while the Python axial element projects that current
onto its rotated axis. The coupled moving-current trajectory is not paired;
the measured differences are in [PARITY.md](PARITY.md).

For pure heave, define only the heave coordinate. Still-water free decay
uses `NoWave()` and `radiation_memory=15`. These moving-body layouts have
no PTO; regular waves require zero heading and a body centered at x=y=0.
The three-DOF Morison force uses a proper pitch rotation and does not copy
the pinned MATLAB source function's nonorthogonal rotation. Other moving
Morison layouts, moving-current wave types, and the normal/tangential coefficient
mode remain unsupported.

For a heave or other `linear_subspace` device, `JONSWAPWave(2.5, 8,
seed=1, gamma=3.3)` selects a JONSWAP sea. Omitting `gamma` uses WEC-Sim's
height/period-dependent value. The published Sphere MPC case uses 2.5 m and
8 s, which infers `gamma=1`; its JONSWAP spectrum therefore equals the PM
spectrum for those inputs. `phase_file` replays a saved MATLAB realization.
For a hydrodynamic body, `PMWave` and `JONSWAPWave` also accept
`frequency_range=(0.5, 1.5)` in rad/s to narrow the BEM frequency interval.
Limits outside the HDF5 range are replaced by that range's endpoints, as in
the pinned MATLAB wave class. A `water_depth` override is currently available
only for fixed Morison bodies; hydrodynamic bodies use their HDF5 depth.
The published Sphere MPC controller has a focused Python runner:

```python
from wecsim import run_sphere_mpc

result = run_sphere_mpc(
    "sphere.h5", "coeff.mat", seed=1,
    max_force=2e6, max_force_rate=1.5e6,
)
heave = result.position
pto_force = result.pto_force
```

This reproduces the published single-body heave layout, fourth-order
radiation fit, 0.5 s optimizer updates, and 0.5 s command transition. The
force, force-rate, heave, speed, horizon, and penalty settings are adjustable.
An integer seed produces a reproducible Python sea; pass a 500-by-1 `phase`
array to replay MATLAB's sea. The paired gate covers the published defaults,
not arbitrary changes to the controller or a general multi-body MPC.

For the published Sphere `Mean_Drift` application, select the control-surface
coefficient in the generated HDF5 file and use convolution radiation:

```python
from wecsim import RegularCICWave, WEC

wec = WEC("Sphere mean drift")
sphere = wec.body("sphere", "Mean_Drift/hydroData/sphere.h5",
                  inertia=(837.75804096,) * 3,
                  mean_drift="control_surface")
for axis in ("surge", "heave", "pitch"):
    wec.coordinate(axis, sphere.move(axis))
result = wec.run(RegularCICWave(0.1, 2), dt=0.01, end_time=100,
                 ramp_time=20, radiation_memory=10)
drift_force = dict(result.raw.extra_outputs)["body1_mean_drift_force"]
excitation_force = dict(result.raw.extra_outputs)["body1_excitation_force"]
```

`mean_drift="momentum_conservation"` selects that HDF5 coefficient when
present. An absent selected dataset raises an error. The regular-wave drift
force is the selected coefficient times wave amplitude squared and the wave
force ramp. The published unmoored linear Sphere run travels more than 16 m
in surge over 100 s; the paired result reproduces that source model but does
not establish physical accuracy at such a large displacement.

The runner rejects unsupported layouts and settings. For PM or JONSWAP, supply
`wave.height`, `wave.period`, optional `directions` and `spreading`, and either
an integer `seed` or a `phase_file` CSV to replay a MATLAB realization.
JONSWAP also accepts an optional positive `gamma`.
For `linear_subspace`, `wave.frequency_range` optionally narrows the HDF5
frequency interval in rad/s. A `wave.water_depth` override is supported only
for fixed Morison bodies.
`simulation` accepts `dt`, `end_time`, optional `ramp_time`, `rho`, `g`, and
`radiation_memory` for radiation-memory cases. For an RM3 `regularCIC` floating
joint, set `simulation.radiation_method` to `"fir"` to use the published
discrete FIR calculation; `"convolution"` remains the default. The same
setting is available for no-wave free decay, without a paired FIR baseline. Both use
`radiation_memory` (60 s by default). The general dynamics module
assembles the supported body and PTO forces; it does not parse Simscape models.
For a `linear_subspace` device, the Python API also accepts a WEC-Sim
three-column imported spectrum MAT file, including its saved phases:

```python
from wecsim import ImportedSpectrumWave

result = wec.run(ImportedSpectrumWave("spectrumData1.mat"),
                 dt=0.1, end_time=400, ramp_time=100,
                 radiation_memory=60, base_dir="path/to/inputs")
```

The file path resolves from `base_dir`. This selects an incident sea for the
configured Python WEC; the specialized RM3 floating-joint MCR runner remains
the paired solver for the published four-coordinate RM3 sea-state motion.
For a sampled time/elevation MAT record, use the same builder with
`ImportedElevationWave`:

```python
from wecsim import ImportedElevationWave

result = wec.run(ImportedElevationWave("etaData.mat"),
                 dt=0.01, end_time=40, ramp_time=10,
                 radiation_memory=15, base_dir="path/to/inputs")
```

The named MAT variable defaults to `etaData` and must contain increasing time
and elevation columns. `reapply_force_ramp=True` explicitly reproduces the
second force ramp used by the pinned MATLAB body block in paired comparisons.
The inherited `WaveClass` now uses the pinned MATLAB PM and JONSWAP spectrum
definitions, including height-dependent PM energy and JONSWAP's inferred
`gamma` when it is unspecified. Its seeded phases use a local NumPy generator
by default. Set `phaseGenerator="matlab"` to reproduce a pinned WEC-Sim
Threefry substream without a saved phase file. A supplied `phaseData` matrix
still replays a specific frequency-by-direction realization:

```python
from wecsim.waveClass import WaveClass

wave = WaveClass("irregular")
wave.T, wave.H = 8, 2.5
wave.spectrumType = "PM"
wave.freqDisc, wave.numFreq = "Traditional", 64
wave.waveDir, wave.waveSpread = [0, 30, 90], [0.1, 0.2, 0.7]
wave.phaseSeed, wave.phaseGenerator = 1, "matlab"
wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)
```

For a regular sea, sample surface elevation at any world XY locations using
the public wave object. The published OSWEC wave-marker case has a 10.9 m BEM
water depth:

```python
import numpy as np
from wecsim import RegularWave

time = np.arange(0, 400.1, 0.1)
locations = [[0, 0], [10, 0], [0, -20]]
markers = RegularWave(2.5, 8).elevation_at(
    time, locations, water_depth=10.9, ramp_time=100,
)
# markers has one row per time and one column per location.
```

The historical BS fixture remains a compatibility check for the
original Python port; current MATLAB WEC-Sim rejects BS inputs.
For this RM3 convolution layout, the solver follows the published pitched
slider geometry and uses an implicit effective added mass by default. To
reproduce the pinned MATLAB/Simulink numerical trajectory, set
`simulation.added_mass_scheme` to `"simulink_delay"` (or pass that keyword to
`solve_rm3_regular` and the RM3 MCR runners). This explicitly selects the
source model's mass split and 1e-7 s acceleration delay; the delay is a
numerical setting, not a WEC property. The published Cases 5 and 6 use a
suspect fitted state-space radiation model and remain unsupported as coupled
trajectories. The paired Cases 3 and 4 tests check the ordinary implicit solver
separately from this opt-in numerical comparison; neither path uses the Cases
5 and 6 fit.

For a prescribed two-body velocity history, the diagnostic
`wecsim.radiation.replay_rm3_fitted_radiation` reproduces the pinned MATLAB
fit's force convention. It does not couple that nonpassive fit into a WEC run.

For the `fixed_hinge` layout, an optional second body can be fixed. The
regular-wave case accepts a nonhydrodynamic base declared with `nonhydro:
true`, `fixed: true`, and a three-component `center_gravity`; PM cases can
instead use a hydrodynamic base with `fixed: true`. Its stationary motion
appears in the response. With this base, `constraint.location` is its ground
attachment and `pto.location` is the flap's hinge attachment. The base's
constraint reaction forces are not yet calculated.
The `regularCIC` floating-joint path has paired MATLAB checks for RM3
body-to-body Cases 3 and 4 and all eight physical settings in the published
RM3 Multiple Condition Runs Option 1 sweep. The published RM3 radiation
options case also checks constant, convolution, and FIR dynamics over 500 s.
MATLAB's fitted radiation state-space Cases 5 and 6 remain outside that
validated path.
The programmatic RM3 solver also accepts `mooring_surge_stiffness` for a
linear spring at its floating joint and `excitation_force` for sampled imported
waves. The pinned `MooringMatrix` case pairs these forces and the full 400 s
body, PTO, and mooring trajectories with MATLAB. Its source body block applies
the wave ramp again after the imported-elevation convolution; the paired test
reproduces that source-specific step explicitly. The same setup can be passed
to the public Python case runner without precomputing a force array:

```python
from wecsim import run_case

case = {
    "simulation": {"dt": 0.01, "end_time": 400, "ramp_time": 40,
                   "radiation_memory": 60},
    "wave": {"type": "elevationImport", "file": "Mooring/MooringMatrix/etaData.mat",
             "variable": "etaData", "reapply_force_ramp": True},
    "bodies": [
        {"hydro_file": "_Common_Input_Files/RM3/hydroData/rm3.h5",
         "hydro_body": 1, "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
        {"hydro_file": "_Common_Input_Files/RM3/hydroData/rm3.h5",
         "hydro_body": 2, "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
    ],
    "constraint": {"kind": "floating_joint", "location": [0, 0, 0],
                   "initial_coordinate": {"spar_heave": -0.21}},
    "pto": {"kind": "relative_heave", "damping": 1_200_000},
    "mooring": {"kind": "joint_surge_spring", "stiffness": 100_000},
}
response = run_case(case, base_dir="path/to/WEC-Sim_Applications")
```

`reapply_force_ramp=True` reproduces the pinned MATLAB body block's second
force ramp; its default is `False`. The current mooring setting is a surge
spring at this joint. Arbitrary mooring matrices and attachment locations
remain unsupported.

For a native MoorDyn model, install the matching
[WEC-Sim MoorDyn library](https://github.com/WEC-Sim/MoorDyn) separately and
pass its library and line-input paths to Python:

```python
import numpy as np
from wecsim import MoorDyn

with MoorDyn("path/to/libmoordyn.so", "path/to/lines.txt").start(
        np.zeros(6), np.zeros(6)) as mooring:
    force_and_moment = mooring.step(body_pose, body_velocity, time, dt)
```

The pose and velocity are six-component vectors at the MoorDyn connection.
For the RM3 floating joint, attach that connection to the spar in body-local
coordinates and run the coupled device through the Python interface:

```python
from wecsim import ImportedElevationWave, MoorDyn, WEC

wec = WEC("RM3 MoorDyn")
float_body = wec.body("float", "rm3.h5", inertia=(0, 21_306_090.66, 0))
spar = wec.body("spar", "rm3.h5", inertia=(0, 94_407_091.24, 0))
wec.floating_joint(
    float_body, spar, damping=1_200_000,
    moordyn=MoorDyn("path/to/libmoordyn.so", "path/to/lines.txt"),
    moordyn_point=spar.at(0, 0, 21.5),
)
result = wec.run(
    ImportedElevationWave("path/to/etaData.mat", reapply_force_ramp=True),
    dt=0.01, end_time=400, ramp_time=40, radiation_memory=60,
    initial_coordinate={"spar_heave": -0.21},
)
connection_force = dict(result.raw.extra_outputs)["moordyn_connection_force"]
```

The solver predicts the connection pose, advances native MoorDyn once per
time step, then applies its force to the two-body motion solve. Its legacy
interface allows one active native model per process. The pinned 400 s RM3
MoorDyn comparison checks both bodies, the PTO, connection force, and three
line tensions. This coupling is currently limited to the reduced RM3 floating
joint with convolution radiation.

The RM3 floating-joint solver accepts optional PTO hard stops in Python:

```python
from wecsim import LinearHardStops, solve_rm3_regular

stops = LinearHardStops(
    lower_bound=-0.6, upper_bound=0.6,
    lower_stiffness=100_000_000, upper_stiffness=100_000_000,
)
response = solve_rm3_regular(
    "rm3.h5", pto_hard_stops=stops, dt=0.025, end_time=400,
)
print(response.pto_stroke, response.pto_stop_force)
```

The case runner also accepts these names under `pto.hard_stops` for a
`floating_joint` with regular waves. Hard stops select adaptive integration,
constant-frequency radiation, and implicit added mass. Ordinary RM3 cases
retain their existing solver. The published MATLAB End_Stops run is
time-step sensitive after contact. Paired motion, force, and energy checks
pass through 400 s against the 0.0125 s MATLAB run; a 0.025 s run gates source
time-step convergence. See
`PARITY.md` for the limits and the published 0.1 s discrepancy.

RM3 multiple-condition runs can be configured in Python without a JSON input:

```python
from wecsim import mcr_grid, run_rm3_mcr

conditions = mcr_grid(
    heights=[1.5, 2.5], periods=[6, 8],
    damping_values=[1_200_000, 2_400_000],
)
result = run_rm3_mcr("rm3.h5", conditions)
power = result.power_matrix(damping=1_200_000)
print(power.periods, power.heights, power.absorbed_power)
first_trajectory = result.traces[0].response
```

`mcr_wave_statistics(path, damping_values)` reads the published Option 2
Excel grid; `mcr_mat_file(path)` reads the Option 3 MAT-file case table.
`run_mcr(conditions, simulate, averaging_start_time=...)` accepts a Python
callback for another configured WEC. The callback returns an `MCRTrace` with
time samples and signed absorbed PTO power; an optional response object keeps
the full trajectory available. The published RM3 example averages from
199.9 s through 400 s and reports positive absorbed power, the opposite of
MATLAB's signed PTO power column. A spring may return stored energy, so an
individual absorbed-power sample may be negative.

The paired MATLAB workflow executes the actual Option 1, 2, and 3
`wecSimMCR` drivers separately. For each option it compares all eight body
and PTO trajectories, average powers, and power matrices. The Option 1
physical conditions are also paired as scalar MATLAB runs. Phase-seed sweeps,
multiple PTOs, and other MCR postprocessing remain unverified.

For the published three imported-spectrum RM3 sea states, pass the MAT-file
table to the Python runner. Paired MATLAB checks cover the incident waves,
excitation, body and PTO trajectories, and absorbed power. Across the three
400 s runs, the largest position differences are 14.4 mm surge, 1.03 mm
heave, and 0.000638 rad pitch; mean absorbed powers differ by at most 29 W:

```python
from wecsim import run_rm3_spectrum_mcr

sea_states = run_rm3_spectrum_mcr(
    "rm3.h5", "RM3_MCROPT3_SeaState/mcrExample.mat",
)
print(sea_states.mean_absorbed_power)
first_wave = sea_states.wave_elevation[0]
first_float = sea_states.traces[0].response.body_position[:, 0, :]
```

The three `spectrumData*.mat` files must sit beside `mcrExample.mat`. Their
third column supplies the exact phase used for each wave component; no random
seed is needed. The runner uses the same 60 s convolution radiation model as
the validated RM3 cases 3 and 4.

For a no-wave floating joint, `constraint.initial_coordinate` and
`constraint.initial_speed` accept four values or dictionaries keyed by
`surge`, `float_heave`, `spar_heave`, and `pitch`. In the published RM3 PTO
extension examples, `float_heave: 5` or `spar_heave: -5` creates the same
initial +5 m relative PTO stroke while moving a different body.
For `linear_subspace`, the map's rows are surge, sway, heave, roll, pitch,
and yaw; its columns are independent generalized coordinates. A two-body
heave case, for example, maps the first body's heave to coordinate 1 and the
second body's heave to coordinate 2. `constraint.initial_coordinate` and
`initial_speed` set those coordinates, and `pto.damping_matrix` and
`stiffness_matrix` apply generalized linear forces. This layout assumes
small rotations and a constant coordinate map.
The matrices may be signed to represent active linear feedback. For example,
the published Sphere reactive PI controller uses a heave-coordinate PTO with
`stiffness_matrix: [[-573350]]` and `damping_matrix: [[49181]]`, giving
`force = -49181 * heave_speed + 573350 * heave_displacement` in the runner's
force convention. This published case has no PTO stroke limit; the paired
trajectory comparison does not validate hardware feasibility.

### Configuring body motions and PTO attachments

The [configurable two-body example](examples/configurable_rm3_pto.json)
shows a WEC definition with named bodies, independent motion coordinates, and
a PTO connection. Run it with:

```sh
python -m wecsim examples/configurable_rm3_pto.json --output results/configurable-rm3.csv
```

In `constraint.coordinates`, each named coordinate lists the body motions it
drives. Assigning the same coordinate to two bodies gives them a shared motion;
separate coordinates let them move independently. A motion names one of surge,
sway, heave, roll, pitch, or yaw and may include a `scale` (default `1`).
For a rotational motion, `pivot` can specify `{"world": [x, y, z]}` or
`{"point": [x, y, z]}` in that body's local center-of-gravity frame. The
rotation then moves the body center about that point. With no `pivot`, the
body rotates about its center of gravity. The example's shared pitch rotates
both RM3 bodies about world `[0, 0, 0]`.
`initial_coordinate` and `initial_speed` may be vectors or dictionaries keyed
by coordinate name. The original 6-by-N `coordinate_map` remains available for
advanced cases. These maps describe constant, small-motion kinematics, not an
arbitrary multibody joint solver.

`ptos` is a list of linear actuators. Each PTO has a unique `name`, `from` and
`to` endpoints, and damping and/or stiffness. A body endpoint uses
`{"body": "float", "point": [x, y, z]}`: the point is in that body's local
frame relative to its HDF5 center of gravity. A fixed world anchor uses
`{"ground": [x, y, z]}`. `axis` gives a fixed world-space force direction;
if omitted, the axis follows the line between endpoints **at the reference
pose** and then stays fixed. Stroke is projected displacement along that axis
from the reference pose. Moving an attachment point changes its moment arm
when the body rotates. The CSV includes each PTO's stroke, velocity, force,
and power absorbed by its damper, as well as named coordinate motion. This
force is positive along the axis on the `to` endpoint and opposite on `from`.
The reported absorbed power is damping times stroke velocity squared for a
constant damper and zero while a declutching PTO is disengaged. This
connection model is linearized for small rotations; it does not update the
actuator's axis as its endpoints move. Use either `ptos` or the older `pto`
matrix in one case.
In the example, each body-local PTO point has a vertical coordinate that
places it at world `z = 0` in the reference pose; its horizontal offset sets
the pitch moment arm.

### PTO tuning

The `fixed_hinge` and `floating_joint` layouts accept scalar `pto.damping`
and optional `pto.stiffness`. They also accept either
`pto.equilibrium_position` or `pto.pretension`. For hinge pitch, the position
is an angle in radians and forces are torques; for the floating joint, it is
relative heave in meters and forces are newtons. The implemented law is
`F = -damping * velocity - stiffness * (position - equilibrium_position)`.
As in MATLAB WEC-Sim, `pretension` sets the equivalent equilibrium position
to `-pretension / stiffness`; nonzero pretension or equilibrium position
requires positive stiffness. These fields default to zero, preserving the
paired reference cases. For example, add
`"stiffness": 100000, "equilibrium_position": 0.1` to the RM3 example's
`pto` object to shift its neutral relative heave by 0.1 m.

The `linear_subspace` layout uses `pto.stiffness_matrix` and
`pto.damping_matrix`, with optional `pto.equilibrium_coordinate` (N values).
Its force law is `F = -K @ (q - q_eq) - C @ q_dot`. A nonzero offset must
produce a spring force. The zero-offset cases have paired MATLAB checks;
nonzero offsets have analytical and case-level tests. Each actuator in `ptos`
can instead set `damping`, `stiffness`, and either `equilibrium_position` or
`pretension`. Connection geometry and nonzero offsets have case-level and
analytical checks. A paired MATLAB Sphere case covers nonzero stiffness,
extra damping while specifying a body-local attachment shifted 1 m in x.
Its heave-only motion cannot validate attachment-location dynamics or a
rotational moment arm.
Force limits, general controller networks, and general hydraulic PTO models
are not implemented. The published RM3 and OSWEC hydraulic layouts have the
focused Python runners described below.

### Reactive direct-drive PTO

The published `Controls/ReactiveWithPTO` Sphere case couples PI control to a
simple direct-drive generator. Its winding inductance and resistance set the
generator torque response time; gear ratio reflects drivetrain inertia and
friction into the heave equation. Configure those properties on a PTO:

```python
from wecsim import RegularWave, SimpleDirectDrive, WEC, WorldPoint

wec = WEC("reactive sphere with direct drive")
sphere = wec.body("sphere", "sphere.h5", mass="equilibrium")
wec.coordinate("heave", sphere.move("heave"))
wec.pto(
    "PTO1", WorldPoint(0, 0, -2), sphere.at(0, 0, 0), axis=(0, 0, 1),
    direct_drive=SimpleDirectDrive(
        kp=380_000, ki=-152_000, torque_constant=7.186,
        gear_ratio=100, drivetrain_inertia=2, drivetrain_friction=1,
        winding_resistance=.483, winding_inductance=5.223e-3,
    ),
)
result = wec.run(RegularWave(height=2.5, period=9.6664),
                 dt=.01, end_time=200, ramp_time=50)
drive = result.ptos["PTO1"].direct_drive
# drive.current, drive.voltage, drive.resistance_loss,
# drive.electrical_power, drive.shaft_torque, drive.mechanical_power
```

`electrical_power` follows the pinned Simulink block's signed `voltage *
current + resistance_loss` output. `result.ptos["PTO1"].absorbed_power` is
positive when mechanical power enters the PTO. Current validation covers one
pure-heave body, zero-heading regular waves, and one vertical direct-drive
connection. The model does not impose voltage, current, or stroke limits.

### RM3 direct linear generator

The published `PTO-Sim/RM3/RM3_DD_PTO` case uses a three-phase linear
generator between the float and spar. Its flux states, electrical angle,
phase currents and voltages, friction, and generated power are available
through a focused two-heave runner:

```python
from wecsim import DirectLinearGenerator, run_rm3_direct_linear_generator

generator = DirectLinearGenerator(
    stator_resistance=4.58, friction=-100, pole_pitch=.072,
    magnet_flux=8, inductance=.285, load_resistance=-117.6471,
)
result = run_rm3_direct_linear_generator("rm3.h5", generator=generator)
print(result.absorbed_power.mean(), result.electrical_power.mean())
```

The negative load resistance and friction follow the source PTO-Sim block's
sign convention. `absorbed_power` and `electrical_power` are positive when
the device absorbs mechanical power and delivers power to the load. This
runner represents the published vertical two-body layout.

The same generator can be configured on body-local PTO endpoints in a
`WEC` with named motions:

```python
from wecsim import DirectLinearGenerator, RegularWave, WEC

generator = DirectLinearGenerator(
    stator_resistance=4.58, friction=-100, pole_pitch=.072,
    magnet_flux=8, inductance=.285, load_resistance=-117.6471,
)
wec = WEC("RM3 with linear generator")
float_body = wec.body("float", "rm3.h5")
spar = wec.body("spar", "rm3.h5")
wec.coordinate("float_heave", float_body.move("heave"))
wec.coordinate("spar_heave", spar.move("heave"))
wec.pto("PTO1", spar.at(0, 0, 0), float_body.at(0, 0, 0),
        axis=(0, 0, 1), linear_generator=generator)
result = wec.run(RegularWave(height=2.5, period=8),
                 dt=.0005, end_time=400, ramp_time=100)
phase_current = result.ptos["PTO1"].linear_generator.phase_current
```

This public configuration is currently limited to constant-radiation
regular waves and cannot be combined with sampled PTO controls. The
attachment points use the existing small-motion, fixed-axis PTO geometry.

### RM3 rectified hydraulic PTO

The published `PTO-Sim/RM3/RM3_cHydraulic_PTO` case has a dedicated Python
runner for its two heave-constrained bodies. Its cylinder, valve, two
accumulators, motor, generator, and PI load are configured as Python objects.
Run `python -m examples.rm3_hydraulic` for a complete example, or use its
`build_pto()` function as a starting point:

```python
from examples.rm3_hydraulic import HYDRO, build_pto
from wecsim import run_rm3_rectified_hydraulic
from wecsim.irregularWave import pm_equal_energy_components

sea = pm_equal_energy_components(
    HYDRO, significant_height=2.5, peak_period=8,
    directions=(0,), spreading=(1,), seed=1,
)
result = run_rm3_rectified_hydraulic(HYDRO, sea, build_pto())
print(result.pto_force.max(), result.generator_voltage.max())
```

The 400 s paired MATLAB test supplies the saved source phases instead of a
Python seed. It checks independent body, hydraulic, and electrical
trajectories. This runner covers the published vertical slider layout;
general hydraulic PTO networks are not yet part of the `WEC` builder. The
published PI law eventually commands negative load resistance, so this
source-specific setting should not be interpreted as a passive load design.

### Instantaneous nonlinear hydrodynamics for heave

The Python builder supports the published `Nonlinear_Hydro/ode4/Regular` and
`ode4/RegularCIC` ellipsoid cases through a heave-only mesh mode. Supply the
BEM HDF5 file and an STL
whose triangle coordinates are relative to the body's center of gravity:

```python
import numpy as np
from wecsim import RegularWave, WEC, WorldPoint

wec = WEC("ellipsoid")
body = wec.body(
    "ellipsoid", "hydroData/ellipsoid.h5",
    geometry_file="geometry/elipsoid.stl",
    nonlinear_hydro="instantaneous", mass="equilibrium",
    inertia=(1_375_264, 1_375_264, 1_341_721),
    drag_coefficient=1, drag_area=np.pi * 5**2,
)
wec.coordinate("heave", body.move("heave"))
wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
        damping=1_200_000)
result = wec.run(RegularWave(height=4, period=6), dt=0.05,
                 end_time=150, ramp_time=50, rho=1025)
```

For convolution radiation, use `RegularCICWave(height=4, period=6)` and
`radiation_memory=60` in `wec.run`.

The mode integrates mesh buoyancy, instantaneous free-surface
Froude–Krylov correction, heave quadratic drag, BEM diffraction/radiation,
and the configured PTO. The STL determines equilibrium mass when
`mass="equilibrium"`; this can differ from the HDF5 displaced volume.
Current validation covers one pure-heave body with its center of gravity at
horizontal origin in zero-direction regular waves, with either constant or
convolution radiation.
The published ode45 variants are tracked separately: MATLAB applies mesh
buoyancy from the preceding 0.05 s sample while the Python mode evaluates it
at the current state. Their motion comparisons are recorded as solver
diagnostics in [PARITY.md](PARITY.md), not as ode45 numerical parity.
Other motions and sea states raise an error until their mesh force and dynamics
checks are paired with MATLAB.

### Variable draft and mass in heave

The published `Variable_Hydro/Variable_Mass` sphere switches among nine BEM
datasets, masses, and draft equilibria every 100 seconds. Define those states
in Python after generating `draft1.h5` through `draft9.h5` with its `bemio.m`:

```python
import numpy as np
from wecsim import HydroState, RegularWave, WEC, WorldPoint

rho = 1025
drafts = range(1, 10)
states = [
    HydroState(
        f"hydroData/draft{draft}.h5",
        mass=rho * np.pi / 3 * draft**2 * (15 - draft),
    )
    for draft in drafts
]
wec = WEC("variable mass sphere")
sphere = wec.variable_body("sphere", states, switch_times=range(100, 900, 100))
wec.coordinate("heave", sphere.move("heave"))
wec.pto("PTO1", WorldPoint(0, 0, 2), sphere.at(0, 0, 0),
        axis=(0, 0, 1), damping=200_000)
result = wec.run(RegularWave(height=1, period=8), dt=0.01,
                 end_time=900, ramp_time=0, rho=rho)
```

Each HDF5 file supplies the corresponding equilibrium center, displaced
volume, added mass, radiation damping, restoring stiffness, and excitation.
The mass and PTO settings come from the Python configuration. Validation covers
one pure-heave body, a regular zero-heading wave, and one passive vertical PTO.
The source model's large late-time excursions are outside the small-motion
range in which linear BEM coefficients can be trusted.

The published Sphere free-decay cases can also be calculated with the focused solver:

```python
from wecsim.linearHeave import solve_heave_free_decay

response = solve_heave_free_decay("path/to/sphere.h5", initial_displacement=1.0)
# response.time, response.position, response.velocity, response.force_total
```

This solver assumes one heave-only body, zero incident waves, and no PTO,
mooring, or nonlinear force. Generate `sphere.h5` with the published
WEC-Sim_Applications Sphere `bemio.m`, or download the HDF5 artifact from the
[MATLAB reference-model run](https://github.com/cmudrc/wec-sim-python/actions/runs/37472556464).
The RM3 regular-wave heave subsystem can be calculated with the same module:

```python
from wecsim.linearHeave import solve_two_body_regular_heave

response = solve_two_body_regular_heave(
    "path/to/rm3.h5", wave_height=2.5, wave_period=8.0,
    pto_damping=1_200_000.0,
)
# response.position and response.velocity each have two body columns;
# response.pto_force is the PTO force acting on body 1.
```

This calculation includes only vertical translation, a linear relative-motion
PTO, and frequency-dependent hydrodynamic coefficients at the incident wave
frequency. The coupled RM3 reference model also predicts both body surge
motions and their shared pitch:

```python
from wecsim.rm3Regular import solve_rm3_regular

response = solve_rm3_regular("path/to/rm3.h5")
# response.body_position and response.body_velocity have shape
# (time_steps, 2 bodies, 6 DOFs); response.pto_force is the heave PTO force.
```

This model uses the published RM3 inertias, regular wave, floating-joint
geometry, and linear PTO. It covers the canonical case's active degrees of
freedom; it does not supply general Simscape joint dynamics.
The OSWEC reference can generate a seeded Python wave realization and pass
its six-component excitation history to the hinge-pitch solver:

```python
from wecsim.hingePitch import solve_hinged_pitch_from_excitation
from wecsim.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)

components = pm_equal_energy_components(
    "path/to/oswec.h5", significant_height=2.5, peak_period=8,
    directions=[0, 30, 90], spreading=[0.1, 0.2, 0.7], seed=7,
)
wave = synthesize_irregular_response(
    "path/to/oswec.h5", components, dt=0.1, end_time=400, ramp_time=100,
)
response = solve_hinged_pitch_from_excitation(
    "path/to/oswec.h5", wave.excitation_force, hinge_z=-8.9,
    body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
)
# response.angle is in radians; wave.elevation is the incident elevation.
```

It models pitch about a fixed hinge and the radiation-memory force. PM
equal-energy bins, directional excitation, and pitch are checked against
current MATLAB WEC-Sim using the same saved random phase matrix. A Python
integer seed creates a reproducible NumPy realization by default; the
`phase_generator="matlab"` option reproduces the pinned source substream.

The two published OSWEC hydraulic PTO applications use a fixed-base flap,
one of two crank linkages, and the configurable rectified hydraulic network.
Run the Python example with a BEMIO OSWEC HDF5 file:

```bash
python -m examples.oswec_hydraulic path/to/oswec.h5
python -m examples.oswec_hydraulic path/to/oswec.h5 --layout fixed
```

`AdjustableRodCrank`, `FixedRodCrank`, `RectifiedHydraulicPTO`, and
`run_oswec_rectified_hydraulic` are also importable from `wecsim`. The example
shows the cylinder, valve, accumulators, motor, generator, and controller
settings as editable Python objects. Paired tests use MATLAB's saved PM
phases; the example uses a reproducible Python seed. This runner covers the
two published pitch-hinge layouts, not arbitrary WEC or PTO attachment
geometry. The adjustable-rod source controller commands negative electrical
load resistance for much of its run, so trajectory agreement does not
establish a passive load design.

The published OSWEC Desalination rod uses a different attachment geometry.
`PitchRodLinkage(anchor=(5.6021271782, -8.7), hinge=(0, -8.9),
body_center=(0, -3.9), body_point=(0.9, -3.1))` specifies its fixed world
anchor and flap-local endpoint in the x/z plane. Its stroke and speed pair
with the 300 s MATLAB source trace.
`ReverseOsmosisMembrane` also exposes the published osmotic valve and linear
hydraulic-resistance settings as a Python object. Given an inlet pressure,
`permeate_flow(pressure)` solves their series pressure drop.
`GasChargedAccumulator` accepts optional atmospheric pressure, dead gas
volume, and hard-stop settings for the published Desalination accumulator.
Given its measured inlet flow and startup pressure, the Python component
integrates liquid volume and reproduces the source pressure trace.
`IdealDoubleActingCylinder(area_a=.26, area_b=.26)` computes the published
incompressible cylinder's rod force and A-in/B-out port flows from chamber
pressures and rod speed. The source force sensor reverses the rod-force sign.
`ReverseOsmosisHydraulicNetwork` predicts high-side pressure and branch flows
from prescribed rod speed, with the source's relief-valve opening lag and
brine-driven motor/pump ratio. Configure it with the published components:

```python
from wecsim import (
    DynamicPressureReliefValve, FourValveRectifiedCylinder,
    GasChargedAccumulator, IdealDoubleActingCylinder,
    ReverseOsmosisHydraulicNetwork, ReverseOsmosisMembrane,
)

network = ReverseOsmosisHydraulicNetwork(
    accumulator=GasChargedAccumulator(
        4, 3e6, atmospheric_pressure=101_325, dead_gas_volume=4e-5,
        hard_stop_stiffness=1e7, hard_stop_damping=1e7,
    ),
    membrane=ReverseOsmosisMembrane(
        .6023e8, 3e6, 5e4, .3, 1e-12, .7, 850,
    ),
    relief=DynamicPressureReliefValve(
        5.6e6, 1e5, .0267, 1e-12, .7, 850, .1,
    ),
    cylinder_area=.26,
    pump_to_motor_displacement_ratio=.95,
    pump_outlet_resistance=.6023e8 / 22,
)
state = network.initial_state()
state = network.step(rod_speed=0.1, previous=state, dt=.01)
print(state.pressure, state.permeate_flow, state.relief_flow)

valves = FourValveRectifiedCylinder(
    cylinder=IdealDoubleActingCylinder(.26, .26),
    max_area=.05, leakage_area=1e-8,
    discharge_coefficient=.7, fluid_density=850,
)
pressure_a, pressure_b = valves.chamber_pressures(0.1, state.pressure)
rod_force = valves.cylinder.force(pressure_a, pressure_b)
```

For the coupled published flap and desalination PTO, pass those configured
`network` and `valves` objects to the runner with the body-local rod endpoint
and Morison elements:

```python
from wecsim import MorisonElement, PitchRodLinkage, run_oswec_desalination
from wecsim.irregularWave import pm_equal_energy_components

hydro_file = "path/to/oswec.h5"
sea = pm_equal_energy_components(
    hydro_file, significant_height=2.64, peak_period=9.86,
    directions=(0,), spreading=(1,), count=250, seed=7,
)
linkage = PitchRodLinkage(
    anchor=(5.6021271782, -8.7), hinge=(0, -8.9),
    body_center=(0, -3.9), body_point=(.9, -3.1),
)
morison = [MorisonElement(
    point=(0, 0, z), drag_coefficient=(1, 1, 1),
    added_mass_coefficient=(0, 0, 0), area=(32.4, 0, 32.4), volume=0,
) for z in (-3, -1.2, .6, 2.4, 4.2)]
response = run_oswec_desalination(
    hydro_file, sea, network, valves, linkage, morison,
)
print(response.pitch, response.high_pressure, response.permeate_flow)
```

The runner advances incident waves, flap motion, radiation memory, Morison
drag, chamber force, and hydraulic pressure together. Its Morison drag follows
the published source convention for this case. For a paired MATLAB comparison,
pass the source's saved phase matrix instead of a Python seed. The raw chamber
force in the pinned Simscape run alternates at the 0.01 s output interval, so
the paired force check uses a 0.1 s mean and a separate total-work check.

The published MBARI Cable application has an axial spring and damper between
two attachment points. Its source force law is available as a Python object:

```python
from wecsim import PlanarCableAttachment, WecSimCableTension, run_mbari_cable

cable = WecSimCableTension(
    stiffness=1_000_000, damping=100,
    length=17.8, initial_length=18,
)
force_z = cable.force_z(relative_position=0.03, relative_velocity=0.1)

attachment = PlanarCableAttachment(
    base_offset=(0, 1.95), follower_offset=(0, -5.2),
    initial_length=18,
)
# Each pose/rate is (x, z, pitch)/(vx, vz, pitch_rate).
stroke, speed = attachment.motion(
    base_pose=(0, -29.95, 0), base_rate=(0, 0, 0),
    follower_pose=(0, -4.8, 0), follower_rate=(0, 0, 0),
)
force_z = cable.force_z(stroke, speed)

# The published application generates this file with Cable/hydroData/bemio.m.
response = run_mbari_cable("applications/Cable/hydroData/mbari.h5",
                           cable=cable, attachment=attachment)
print(response.body_position.shape)  # time × three bodies × six coordinates
print(response.cable_force_z)         # signed axial force, N
```

`relative_position` and `relative_velocity` are the source cable block's
local z displacement and speed. `PlanarCableAttachment` computes these from
body-local endpoints as the bodies surge, heave, and pitch. The coupled runner
advances the two hydrodynamic bodies and the nonhydrodynamic cylinder, including
spring/damper feedback and quadratic drag on the cable's two endpoint bodies.
Its planar buoy-to-cylinder spherical joint matches the published geometry.
It uses implicit regular-wave added mass and omits the endpoint bodies' 1 kg
inertias. The MATLAB source shifts added mass into its Simscape bodies and
feeds the remainder through a short acceleration delay, so brief cable snap
events differ more than body positions. The source law can report positive
force while stretched and contracting rapidly because its damping term is not
clamped.

The published WaveBot `CalcImpedance` application uses three independent
motions and a sampled multisine input. Generate its HDF5 file with
`Load_Mitigating_Controls/hydroData/bemio.m`, then run:

```python
from wecsim import run_wavebot_impedance

response = run_wavebot_impedance(
    "applications/Load_Mitigating_Controls/hydroData/waveBotBuoy.h5",
    "applications/Load_Mitigating_Controls/CalcImpedance/multisine3DOFA.mat",
)
print(response.body_position[:, [0, 2, 4]])  # surge, heave, pitch
```

The default uses positive diagonal damping in those coordinates. MATLAB's
published `linearDamping(1:2:5)` input instead fills the first matrix column,
making all three damping forces depend on surge speed. Use
`source_linear_damping=True` to compare the published source trajectory. The
runner applies the source model's negative pitch-actuation sign in both modes.
The adaptive `ControlTests` application has a separate published-case runner
and paired gate in [PARITY.md](PARITY.md).

The published WaveStar base application can also be run without MATLAB after
generating `WECCCOMP/hydroData/wavestar.h5` with its `bemio.m`:

```python
import numpy as np
from wecsim import run_wavestar_published
from wecsim.irregularWave import jonswap_equal_energy_components

hydro = "applications/WECCCOMP/hydroData/wavestar.h5"
sea = jonswap_equal_energy_components(
    hydro, significant_height=0.0625, peak_period=1.412,
    directions=np.array([0.]), spreading=np.array([1.]),
    gamma=1, seed=1,
)
response = run_wavestar_published(hydro, sea)
print(response.angle, response.pto_stroke)
```

By default, the function integrates the published unforced, single-coordinate linkage and
uses the source application's fitted state-space radiation model. A Python
seed gives a reproducible sea; replaying MATLAB's exact sea requires its saved
phase column. The fit is positive in the active joint coordinate at sampled
frequencies but differs from the BEM damping table. The broader `WEC` runner's
physical radiation default is unchanged. The fault application has a separate
paired trajectory below; nonlinear predictive closed-loop trajectory parity remains open.

The published nonlinear predictive application also has a paired PTO actuator
component. It accepts the controller's sampled torque request and the current
rod stroke every 0.05 s:

```python
from wecsim import (
    WaveStarNmpcActuator, WaveStarNmpcController,
    WaveStarNmpcObserver, WaveStarNmpcPredictor,
)

actuator = WaveStarNmpcActuator()
axial_force = actuator.step(command_torque, stroke)

observer = WaveStarNmpcObserver(dt=0.05)
predictor = WaveStarNmpcPredictor()
controller = WaveStarNmpcController()
for time, stroke, previous_command_torque in sampled_inputs:
    estimated_state = observer.step(stroke, previous_command_torque)
    future_excitation_moment = predictor.step(estimated_state[-1])
    command_torque = controller.step(
        estimated_state, future_excitation_moment, time,
    )
```

The controller converts torque to force with its own idealized neutral rod
length, which differs slightly from the physical B-to-C linkage. The source
sea, linkage, resistive startup, and actuator force have paired gates. A
derived 0.001 s MATLAB run also pairs the unforced plant through 9.95 s; the
published 0.05 s run has a larger numerical difference. The observer replays
all five published estimator states from saved stroke and the previous torque
request. The autoregressive forecast is paired from saved estimated moments;
the composed observer and predictor are also paired after 11 s from stroke
and the preceding torque request. Its 10–11 s startup fit is numerically
unstable in the source. Given saved estimator states and forecasts, the Python
controller replays every published torque command within `1e-4` N m over the
full 225 s case. Composing the Python observer, predictor, and controller from
saved stroke and preceding torque request meets the same gate. This command
replay uses the saved MATLAB motion as input; independent closed-loop NMPC
motion remains unpaired.

A separate 1 ms MATLAB `ode4` diagnostic runs the same WaveStar sea and
geometry through 14.95 s with resistive control active from 10 s. An
independent Python run pairs pitch within `5e-5` rad, PTO stroke within
`10` µm, and axial PTO force within `0.03` N across all 14,951 samples.
That derived run also samples the controller at 1 ms; the published 0.05 s
`ode8` closed loop and NMPC motion after 15 s remain open.

To run the complete 225 s Python controller from its own motion, pass the
sampled PTO to the WaveStar runner:

```python
import numpy as np
from wecsim import WaveStarNmpcPTO
from wecsim.irregularWave import jonswap_equal_energy_components
from wecsim.wavestar import run_wavestar_published

hydro = "applications/WECCCOMP/hydroData/wavestar.h5"
sea = jonswap_equal_energy_components(
    hydro, significant_height=0.1042, peak_period=1.836,
    directions=np.array([0.]), spreading=np.array([1.]),
    gamma=3.3, seed=1,
)
pto = WaveStarNmpcPTO(plant_dt=0.001, control_dt=0.05)
motion = run_wavestar_published(
    hydro, sea, dt=0.001, end_time=225, ramp_time=25,
    g=9.80665, pto_controller=pto, output_stride=50,
)
print(motion.angle, pto.command_torque, motion.pto_force)
```

Python's seed generates a new realization; replaying the published MATLAB
sea requires its saved phase column. The published 0.05 s `ode8` case and
this fine-step Python run show a bounded but material motion difference, so
their full 225 s trajectory is recorded as a diagnostic rather than a
paired parity gate.

The WaveStar fault application has separately paired PTO components:

```python
from wecsim import StribeckFriction, WaveStarFaultController

joint = StribeckFriction(.25, .1, .2, .001)
joint_torque = joint.torque(angular_speed)
controller = WaveStarFaultController(gain=10, filter_frequency=30)
axial_force = controller.step(stroke, noise=position_noise, dropout=False)
```

Call `step` every 0.001 s with the current B-to-C stroke. Its sensor noise
and dropout inputs make the stochastic fault realization explicit. For the
published fault application's coupled dynamics, use the same sea and HDF5
from above with a chosen sensor disturbance:

```python
from wecsim import run_wavestar_fault_published

rng = np.random.default_rng(1)
samples = 141_201  # 141.2 s at 0.001 s, including the initial sample
noise = rng.normal(0, 0.003, samples)
dropout = rng.random(samples) < 0.03
response = run_wavestar_fault_published(hydro, sea, noise, dropout)
print(response.angle, response.pto_force)
```

For an exact paired MATLAB run, supply its realized noise and dropout record.
The runner calculates its own body motion, sensor geometry, controller states,
and forces. The source's fitted radiation model is confined to this case; see
[PARITY.md](PARITY.md) for the full paired gates and limits.

The paired 300 s test drives the network and valve model with MATLAB's saved
rod speed, without saved pressure, flow, or force as inputs. The legacy source
chamber pressure flips on alternating 0.01 s samples while valve flow stays
smooth. The test therefore compares the resolved force over 0.1 s and total
rod work, while reporting the raw discrepancy. Coupled flap motion remains an
open gap; see [PARITY.md](PARITY.md).

For a supported fixed-hinge OSWEC, use a Python wave object with a
frequency-resolved directional MAT spectrum:

```python
from wecsim import FullDirectionalSpectrumWave, WEC

wec = WEC("OSWEC")
flap = wec.body("flap", "path/to/oswec.h5", mass=127000,
                inertia=(0, 1.85e6, 0), hydro_body=1)
base = wec.body("base", "path/to/oswec.h5", mass=999,
                inertia=(999, 999, 999), hydro_body=2)
wec.fixed_hinge(flap, base, damping=12000)
result = wec.run(
    FullDirectionalSpectrumWave("path/to/fullDirSpectrum.mat", seed=7),
    dt=0.05, end_time=400, ramp_time=100, radiation_memory=30,
)
```

The wave also accepts a saved `phase_file` CSV instead of a seed. Its default
includes each heading-bin width in both elevation and force, preserving the
spectrum's integrated energy. The published MATLAB `Full_Directional_Waves`
force block omits that width even though its elevation includes it. To
reproduce that specific source trajectory, set
`force_quadrature="matlab_omitted"` and
`excitation_interpolation="spline_frequency"` on the wave, and set
`phase_generator="matlab"` with the source seed. The optional
`wec.fixed_hinge(..., added_mass_scheme="simulink_delay")` compares against
Simulink's delayed added-mass feedback. Integrated forcing and implicit added
mass remain the defaults. The equivalent low-level case wave type is
`spectrumImportFullDir`.

The same wave calculation is available directly in Python:

```python
from wecsim import (imported_full_directional_components,
                    synthesize_full_directional_response)

components = imported_full_directional_components(
    "path/to/oswec.h5", "path/to/fullDirSpectrum.mat", seed=7,
)
incident = synthesize_full_directional_response(
    "path/to/oswec.h5", components,
    dt=0.05, end_time=400, ramp_time=100,
)
# incident.elevation and incident.excitation_force are NumPy arrays.
```

For a regular-wave ParaView surface, the legacy wave object can write numbered
VTP wave meshes and `ground.txt` using WEC-Sim's grid layout:

```python
from wecsim.waveClass import WaveClass
from wecsim.paraviewClass import ParaviewClass

waves = WaveClass("regular")
waves.H, waves.T, waves.waveDir = 2.5, 8, [0]
waves.waveSetup([0.3, 1.5], 30, 0, 0.1, 21, 9.81, 1000, 2)
ParaviewClass(waves).write_paraview_vtp_wave(
    [0, 1, 2], "results/paraview", domain_size=40,
)
```

For a body mesh, `ParaviewClass(waves).write_paraview_vtp(...)` accepts
zero-based triangular `faces`, local `vertices`, one six-DOF `pose` per time,
and optional face-pressure arrays. It writes numbered body VTP files with the
source XYZ rotation order. The published OSWEC and RM3 visualization cases
still need direct output gates for their actual meshes and mooring lines.

To run a supported case without writing Python code:

```sh
python -m wecsim.reference rm3 --h5 path/to/rm3.h5 --output results/rm3.csv
python -m wecsim.reference rm3 --h5 path/to/rm3.h5 --output results/rm3-b2b.csv --b2b
python -m wecsim.reference oswec --h5 path/to/oswec.h5 --output results/oswec.csv --seed 7
python -m wecsim.reference sphere --h5 path/to/sphere.h5 --output results/sphere.csv --initial-displacement 1
```

Each command writes a numeric CSV and an adjacent JSON file with the model
settings, HDF5 SHA-256 hash, NumPy version, and Git revision/dirty state.
For RM3, `--b2b` includes cross-body hydrodynamic coupling as in the
published B2B Case 2; the default matches B2B Case 1.
The original README projected completion in
August 2022; that date is no longer applicable. See [PARITY.md](PARITY.md)
for the tested scope and next reference case.

## Design direction

The Python dynamics engine uses independent coordinates for each supported
constraint layout, avoiding numerical joint drift. New layouts and force
models need explicit MATLAB comparisons before they are called supported.
Wave and body VTP serialization now have small pinned MATLAB source gates;
full application visualization and BEMIO conversion remain future work.
