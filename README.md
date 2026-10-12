# wec-sim-python

[![Published-case alignment](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2Fcmudrc%2Fwec-sim-python%2Fmaster%2Fdocs%2Fparity-cases.json)](PARITY.md#reference-model-and-case-inventory)
[![Python tests and coverage](https://github.com/cmudrc/wec-sim-python/actions/workflows/python-parity.yml/badge.svg?branch=master)](https://github.com/cmudrc/wec-sim-python/actions/workflows/python-parity.yml)
[![Selected MATLAB checks](https://github.com/cmudrc/wec-sim-python/actions/workflows/reference-model-baselines.yml/badge.svg?branch=master&event=push)](https://github.com/cmudrc/wec-sim-python/actions/workflows/reference-model-baselines.yml)

A Python library for configuring and simulating wave energy converters, with
paired checks against MATLAB [WEC-Sim](https://github.com/WEC-Sim/WEC-Sim).
This project builds on Sungjun Won's original WEC-Sim-Python work; his commits
remain in the repository history.

This is a pre-alpha parity effort. Supported Sphere, RM3, OSWEC, barge, and
selected application cases have MATLAB checks, but the library does not yet
cover all WEC-Sim models or physics. The [parity ledger](PARITY.md) records
which results are independently paired and what remains open.

## Install

Requires Python 3.12 or newer. Install from a clone:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e .
```

Or install from GitHub:

```sh
python -m pip install "git+https://github.com/cmudrc/wec-sim-python.git"
```

Both methods provide the `wecsim` Python package and command.

## Quick start

```python
from wecsim import NoWave, WEC, WorldPoint

wec = WEC("Heaving float")
body = wec.body("float", "path/to/hydro.h5")
wec.coordinate("heave", body.move("heave"))
wec.pto("main", body.at(0, 0, 0), WorldPoint(0, 0, 10), damping=1_200_000)

result = wec.run(NoWave(), dt=0.1, end_time=20,
                 radiation_memory=15, initial_coordinate={"heave": 1})
heave = result.bodies["float"].position[:, 2]
pto_force = result.ptos["main"].force
```

`from wecsim import ...` is the recommended import pattern. The WEC builder,
configuration, and result types can also be imported from `wecsim.interface`.
See the runnable [two-body RM3 example](examples/configurable_rm3_pto.py) for
body-local PTO attachments and shared coordinates.

## Documentation

- [User guide](docs/guide.md): WEC configuration, wave inputs, PTOs, and application examples.
- [Parity ledger](PARITY.md): pinned MATLAB reference cases, evidence, and open gaps.
- [Examples](examples/): runnable cases and bundled inputs.

The source is organized by role under `wecsim/`: `interface`, `environment`,
`hydrodynamics`, `dynamics`, `components`, `control`, `analysis`, and
`visualization`.

## Tests and coverage

```sh
python -m pip install -e ".[test]"
python -m pytest -q --cov=wecsim --cov-branch --cov-report=term-missing
```

[Python CI](https://github.com/cmudrc/wec-sim-python/actions/workflows/python-parity.yml)
runs the suite on each push and pull request, summarizes line and branch
coverage, and uploads `coverage.xml` and `coverage.json`. Many MATLAB-dependent
tests skip without their reference assets. The bundled
[RM3 MATLAB trajectory](tests/fixtures/matlab/README.md) runs in ordinary CI.
Code coverage measures executed Python paths; the [parity ledger](PARITY.md)
measures agreement with WEC-Sim.

## License

See [LICENSE](LICENSE).
