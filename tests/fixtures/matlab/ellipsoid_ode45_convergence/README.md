# Ellipsoid ode45 derived source-step screen

This fixture retains full 150 s R2025b source histories for the
`Nonlinear_Hydro/ode45/Regular` and `RegularCIC` applications. The source
workflow checks out WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and WEC-Sim_Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. It creates a copy of each
application input for every run, replacing the published `simu.dt = 0.05`
with 0.005, 0.0025, 0.00125, or 0.000625 s. At runtime it verifies that
`simu.dt`, `nonlinearDt`, `cicDt`, and Simulink `MaxStep` all equal the chosen
physical step. The first three arms also use that step for `dtOut`. The finest
arm explicitly sets `dtOut=0.005` s and retains 0.000625 s for all four
physical clocks. Each arm runs in a fresh R2025b process over the full 150 s.
The workflow generates the archived BEM HDF5 from pinned WAMIT output in the
application checkout and copies its STL; neither upstream source repository
is changed.

`source.npz` stores all eleven output channels and five checked source settings
at every source step. The finest history uses the common 0.005 s output grid.
`hydroData/ellipsoid.h5` and `geometry/elipsoid.stl` let the portable Python
test advance an independent physical model from the same input. The archive
tool checks all source settings, validates the 0.005 s history against the
earlier [full-step source archive](../ellipsoid_ode45_maxstep_only/README.md),
checks wave identity, and checks contraction of heave, speed, PTO force, PTO
power, and signed PTO work across the refinement sequence. Pointwise source
changes are screened on the same 0.005 s grid, with their original numerical
limits unchanged. Signed work is integrated on the full available grid for
the first three arms and on the 0.005 s grid for the finest; the test checks
the full/common integration difference in each earlier arm.
For every arm, the archive and portable gate also check PTO vertical speed
against body heave speed and verify `Fz = -1,200,000 vz` and `P = -Fz vz`
within 1e-9 m/s, 1e-6 N, and 1e-6 W, respectively. These identities reject
the previously encountered inactive local-axis export.

The [controlled output-grid A/B](https://github.com/cmudrc/wec-sim-python/actions/runs/40104326372)
found exact common-grid trajectory and PTO equality when `dtOut` changed from
0.00125 to 0.005 s with the physical clocks fixed at 0.00125 s. This supports
an output-sampling explanation but establishes neutrality only at that tested
physical step. Four all-clock 0.000625 s attempts, including an isolated fresh
process, lost hosted runners during simulation without a MATLAB or physical
gate result. The mixed-output finest source arm is therefore labeled as a
derived source input in its own right, with unchanged full physical gates.

The source model retains its 1e-7 s added-mass feedback Transport Delay. Its
logged restoring force samples the preceding state at the 0.00125 s all-clock
arm. The 0.000625 s arm does not log every solver step, so its immediately
preceding source state is unavailable. The independent Python model evaluates
force from its current state; the test exposes this timing law explicitly at
0.00125 s. Passing this step-screened, derived reference does not change the
separate, pinned 0.05 s published-input
diagnostic or its `qualified_motion` status.

The portable gate rechecks all source convergence conditions, advances Python
at 0.0025 and 0.00125 s over the full duration, and compares its finest run
directly against the 0.000625 s source on the 0.005 s output grid. Its original
pointwise motion/PTO limits remain unchanged. The Python self-step limits are
0.1 mm, 0.1 mm/s, 100 N, and 100 W. The signed source-work final-step change
must be below 0.1% and contract by a factor below 0.7 at each refinement;
the direct Python/source signed-work difference must also be below 0.1%.

Original all-clock source run: [R2025b workflow 40100441636](https://github.com/cmudrc/wec-sim-python/actions/runs/40100441636).
Isolated all-clock diagnostic: [R2025b workflow 40105605489](https://github.com/cmudrc/wec-sim-python/actions/runs/40105605489).

The archived corrected-axis source histories came from these Python-repository
workflow runs and exact commits:

| Source arms | R2025b run | Python-repository commit |
| --- | --- | --- |
| Both cases, 0.005/0.0025/0.00125 s | [40107175117](https://github.com/cmudrc/wec-sim-python/actions/runs/40107175117) | `daa282c9c3ccb6041b7ccccb89304944ea351d36` |
| Regular, physical 0.000625 s with `dtOut=0.005` s | [40106411780](https://github.com/cmudrc/wec-sim-python/actions/runs/40106411780) | `44b0809b295e2163d433be98b246be2ebf58113d` |
| RegularCIC, physical 0.000625 s with `dtOut=0.005` s | [40107064277](https://github.com/cmudrc/wec-sim-python/actions/runs/40107064277) | `1e13e500f815dcb91168e985c2af060bc95b077d` |

The measured source changes below compare adjacent arms on their common
0.005 s grid. Each channel contracted through both intermediate steps, with
successive change ratios below the predeclared 0.7 bound.

| Case | Final heave change | Final speed change | Final PTO force change | Final PTO power change | Final signed-work change |
| --- | ---: | ---: | ---: | ---: | ---: |
| Regular | 0.000231452 m | 0.000341197 m/s | 409.436 N | 333.393 W | 0.069817% |
| RegularCIC | 0.000241166 m | 0.000362903 m/s | 435.484 N | 360.546 W | 0.073934% |

The four signed-work measurements on each arm's available output grid are
20.154279, 20.098017, 20.069978, and 20.055976 MJ for Regular; and 20.186951, 20.127222,
20.097460, and 20.082612 MJ for RegularCIC. The final source-work change is
below the declared 0.1% limit for both cases. Common-grid work also contracts
and passes the same limit. Source settings, wave identity, and archived 0.005 s
continuity pass for both cases. All 38 BEM HDF5 datasets agree in value across
the source arms, and the STL files have identical SHA-256 hashes.

The independent Python 0.0025→0.00125 s self-step gates pass. Against the
finest source arm, direct Python heave/speed/PTO-force/PTO-power errors are
0.000231346 m, 0.000338943 m/s, 406.732 N, and 328.493 W for Regular;
0.000241018 m, 0.000362768 m/s, 435.321 N, and 359.533 W for RegularCIC.
Their signed-work differences are 0.069727% and 0.073878%, respectively,
below the unchanged 0.1% direct gate. These measurements apply only to this
derived source-step screen.

Archived file SHA-256 hashes:

| File | SHA-256 |
| --- | --- |
| `source.npz` | `8dd1d8d28feab13617eab78afad02e128e1fee017364db906af7d29f4e808f23` |
| `hydroData/ellipsoid.h5` | `cd437b3fcdb2837b1534abadc3591d21cf3f001812407973b1cdd12517c5f489` |
| `geometry/elipsoid.stl` | `a225213457fb9273a327a141273b9d68266e650cf973a98f84e84e8607cbb7ec` |
