# JSON `simulation` validation pilot (#634)

**Decision:** Keep the existing shared `_section` and `_number` validators. Do not
add Pydantic as a runtime dependency for this boundary.

## Bounded prototype

The current `simulation` boundary is nine nonblank lines in
`wecsim/dynamics/caseDynamics.py` (the `_section` call and five `_number`
calls). The two helpers are already used by other case sections, so replacing
this section cannot remove them. An executable local prototype using Pydantic
2.14.1 required 34 nonblank implementation lines for the model and adapter,
before changing the `run_case` call site. That is at least 25 additional lines
for this section
and would add Pydantic plus its transitive runtime packages.

| Measure | Current parser | Pydantic `simulation` prototype |
| --- | ---: | ---: |
| Section validation lines | 9 | 34 plus call-site wiring |
| New runtime packages | 0 | 5 in the pilot environment |
| Shared `_section` / `_number` still needed | Yes | Yes |

The prototype used `extra="forbid"`, `strict=True`, `allow_inf_nan=False`, and
`validate_default=True`, then delegated scalar conversion to `_number` to
preserve numeric strings, boolean rejection, finite ranges, and existing
`ValueError` messages. It called `_section` first to preserve the sorted
missing/unknown-key diagnostic. It kept the optional solver, radiation, and
added-mass settings untyped: these have layout-specific validation and the
`radiation_memory` default varies by layout. The resulting JSON Schema had
nine named fields, but a faithful schema for the accepted numeric-string inputs
would need further customization. Modeling the four physics settings more
strictly would change the input contract unless their layout rules moved too.

The original mapping must remain available to the layout checks. They use
`set(sim)` and key presence; for example, an explicit `ramp_time=100` on a
no-wave free decay is rejected even though omitted `ramp_time` is accepted.
A model dump with defaults would erase that distinction. A future model would
need to retain the mapping or use `model_dump(exclude_unset=True)`.

In a local comparison, the prototype produced identical values,
explicit-key sets, exception types, and messages to the current parser for 13
cases: minimal input, numeric
strings, explicit defaults, booleans, null, nonnumeric input, negative and
nonfinite values, unknown and missing keys, a non-object section, and multiple
invalid fields. The permanent tests in
`tests/test_case_simulation_boundary.py` capture the accepted/rejected contract
and explicit-key behavior without requiring Pydantic.

The rejected prototype was:

```python
from typing import Any
from pydantic import BaseModel, ConfigDict, ValidationError, ValidationInfo, field_validator
from wecsim.dynamics.case_inputs import _number, _section

_REQUIRED = {"dt", "end_time"}
_ALLOWED = _REQUIRED | {
    "ramp_time", "rho", "g", "radiation_memory", "radiation_method",
    "added_mass_scheme", "solver",
}

class SimulationInput(BaseModel):
    model_config = ConfigDict(
        extra="forbid", strict=True, allow_inf_nan=False, validate_default=True,
    )
    dt: float
    end_time: float
    ramp_time: float = 100
    rho: float = 1000
    g: float = 9.81
    radiation_memory: Any = None
    radiation_method: Any = None
    added_mass_scheme: Any = None
    solver: Any = None

    @field_validator("dt", "end_time", "ramp_time", "rho", "g", mode="before")
    @classmethod
    def _numeric(cls, value: Any, info: ValidationInfo) -> float:
        name = info.field_name
        return _number(value, f"simulation.{name}",
                       positive=name in {"dt", "rho", "g"},
                       nonnegative=name in {"end_time", "ramp_time"})

def validate_simulation(value: Any) -> SimulationInput:
    value = _section(value, "simulation", _REQUIRED, _ALLOWED)
    try:
        return SimulationInput.model_validate(value)
    except ValidationError as exc:
        raise exc.errors()[0]["ctx"]["error"] from exc
```

## Verification

On the current Python parser, the case-boundary, CLI/provenance,
builder-recovery, Python API, Sphere fixture parity, and RM3 two-column
fixture parity suites passed without Pydantic installed: **138 tests in
8.38 s**. The command was:

```text
python -m pytest -q tests/test_case_simulation_boundary.py \
  tests/test_case_physical_boundaries.py \
  tests/test_unsupported_case_boundaries.py tests/test_cli_bundled_case.py \
  tests/test_builder_recovery.py tests/test_public_builder_recovery.py \
  tests/test_python_api.py tests/test_case_dynamics.py \
  tests/test_sphere_portable_parity.py tests/test_rm3_two_column_motion.py
```

No Python runtime parser or WEC builder code changed in this pilot.
