# WaveStar NMPC source reference convergence diagnostic

The `source_*.npz` files preserve 4,501 samples of independent 225 s
MATLAB R2025b WaveStar NMPC runs from [Python-repository CI run
39883276159](https://github.com/cmudrc/wec-sim-python/actions/runs/39883276159)
and its [repeat run
39886232300](https://github.com/cmudrc/wec-sim-python/actions/runs/39886232300).
They are **diagnostics, not accepted physical parity references**. Each run's
three source artifacts contain the original body/PTO CSV files and controller MAT
files. `convergence.json` and `convergence_rerun.json` record the full
source-to-source comparisons in the two runs.

The workflow pins WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`,
Applications `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`, and the exact
`wavestar.h5` in
[`../wavestar_nmpc_split_step_225s`](../wavestar_nmpc_split_step_225s/README.md)
(SHA-256 `b2541ff01c319987d9157d1be9fbc9ec39a7660d0572d8e316c428b8aa3359d1`).
The 500 wave components match that pinned fixture exactly. All runs use the
`ode4` source plant, the same 50 ms controller and output interval, and plant
steps of 1, 0.5, or 0.25 ms. Only disposable Applications input files were
configured for the runs; no changes were pushed to either MATLAB repository.
The NPZ arrays are float64 copies of exported source time, pitch, pitch
speed, PTO stroke and speed, axial force, requested torque, and components.

| Comparison | Max pitch gap | Max request gap |
| --- | ---: | ---: |
| Older pinned 1 ms vs fresh 1 ms | 0.010596 rad | 0.792728 N m |
| Fresh 1 ms vs 0.5 ms | 0.009193 rad | 0.573686 N m |
| Fresh 0.5 ms vs 0.25 ms | 0.004763 rad | 0.450018 N m |
| First vs repeated 0.25 ms | 0.012014 rad | 0.732573 N m |

The original pointwise limits are 0.0005 rad pitch and 0.1 N m request.
The older and fresh 1 ms source pitches agree to `2.29e-15` rad through
15.2 s, then separate after NMPC activation. The first 15.2 s pitch gaps
between the fresh 1/0.5 and 0.5/0.25 ms runs are `2.6722e-5` and
`1.3348e-5` rad; smaller plant steps have not delivered an acceptable full
closed-loop target. Independent physical Python, run at the matching plant
step with only the exact HDF5 and generated sea phases driving feedback,
differs from the 0.5/0.25 ms source traces by up to 0.006714/0.009435 rad
pitch. The 0.25 ms pitch-peak aggregate gap is 0.7149%, over its 0.5% limit.
These first matched-step Python checks ran locally on macOS ARM64 with NumPy
2.5.4 and SciPy 1.18.1. In the repeat Linux run, the two 1 ms and two 0.5 ms
source exports reproduced exactly at their respective steps. The 0.25 ms
source exports agree within `1.10e-15` rad pitch through 15.2 s, but differ
by up to `0.012014` rad pitch and `0.732573` N m request over 225 s despite
the same HDF5 hash, sea components, source revisions, and step settings.
The repeated 0.25 ms source is `source_0p25ms_repeat.npz`; its matching Linux
Python pitch gap is `0.011711` rad and the pitch-peak aggregate gap is
`0.7851%`. Both exceed the original limits. The earlier pinned 1 ms trace
also differs from both standalone 1 ms runs, although those standalone runs
match each other numerically. This record does not isolate which runner or
source setup detail caused the difference.

SHA-256:

- `source_1ms_repeat.npz`: `42f3264022177833634731711086cc90a1a39eb313ccdc05110123eb9c3035f6`
- `source_0p5ms.npz`: `71d5745bfd28ad25be4b25fb73b6656d9374c6923d6059d72c1ec6d105e63ad3`
- `source_0p25ms.npz`: `0f62bc5692e845f739bb0df715970d006124d4514cdb928059fb2115263596a6`
- `source_0p25ms_repeat.npz`: `5f6f576b94e4064352ce5bf1ba7a83e69c83faf157c8f1abe77488645b627d24`

The Python repository's
[`wavestar-nmpc-reference-convergence.yml`](../../../../.github/workflows/wavestar-nmpc-reference-convergence.yml)
regenerates the MATLAB exports, audits autonomous physical Python, and
compares the source traces. A repeatable, step-converged source target is
still required before issue #623's pointwise gate can close.
