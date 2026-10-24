# RM3 MoorDyn visualization corrected baseline

This fixture pairs an independent 80 s Python WEC and live native MoorDyn
trajectory with a MATLAB R2025b WEC-Sim source run. The disposable MATLAB
case starts from the pinned published `Paraview_Visualization/RM3_MoorDyn_Viz`
application. It makes exactly two input edits: `simu.dt` changes from 0.01 to
0.0003125 s, and `waves.phaseSeed = 1` fixes the random sea. It retains all
other published settings, including the 0.1 s body output and ParaView
frame grid, `ode45`, Traditional JONSWAP sea, and all MoorDyn settings.
The unchanged 0.01 s published trajectory and its strict expected-failure
test remain separate in `test_rm3_moordyn_viz_source.py`.

## Provenance

- Definitive corrected source run: [GitHub Actions 40049296662](https://github.com/cmudrc/wec-sim-python/actions/runs/40049296662). The earlier [40048188581](https://github.com/cmudrc/wec-sim-python/actions/runs/40048188581) established the corrected 801-frame scene and 0.1 s histories; the definitive run also exports all 8,001 original 0.01 s fairlead samples.
- Pinned WEC-Sim: `0753b2e47f2457c078751dcfe5d251d1767b80ab`.
- Pinned WEC-Sim Applications: `d53d4d4c9eda2581f04204f5d394a6ef84bb099e`.
- Pinned MoorDyn: `34c0703c9a5b74fd2cc5bb52caaadaaa356d675e`.
- Published application input SHA-256: `16a0d1519e218429924bcd879335abdb86b51bbf91fa0f5a0a781d93a6ad9ff5`.
- Disposable corrected application input SHA-256: `f403e48a842fb1c0fe73b73d7476dcf4b652794a449bee08a077f3c66dfe8df1`.
- Corrected source generator SHA-256: `c1d18ccf59793ccfe8d7dffee899f082ac9fd868fe79fc16ba91affb4daa7541`.
- Bundled `rm3.h5` SHA-256: `e0c163a8f2e2254e63fb6261fac8b3e539592efb8f5d4db7ec7e0f5851668070`, copied from the definitive source run. A fresh BEMIO container may have a different byte hash; the hosted gate compares all 84 HDF5 datasets exactly, including body properties, restoring, excitation, and radiation kernels.
- Bundled `lines.txt` SHA-256: `d98ac9edcd8077363fe899f644a88ffc2f12ab88215ac2c35fdf6b6e0af99ded`. The hosted gate also checks byte equality with the pinned Applications input.
- Source body STL SHA-256: `float.stl` `3d990ab8584fa21af2f22278bfe862fe119a491e4caa0d8d9abd988d1096bf4b`; `plate.stl` `9db8cc11ecf74373c0b570651e91e822617379b58f4137432028cfbb1d107acb`. The hosted gate checks these against the repository meshes used by Python.

`source.npz` contains the 0.000625 s and 0.0003125 s source refinement arms and
the corrected 0.0003125 s source scene arm. It retains 801 wave/body/PTO/
connection/line frames, 8,001 fairlead samples for every arm, all line nodes
and segment tensions, and compact signatures parsed from every actual
corrected MATLAB VTP frame. The signatures include wave samples and extrema,
body centroid/bounds/total area, and every line node and segment tension.
Wave and body topology and VTP finite values are checked while compacting.
Its SHA-256 is `faf100346d7016f3c243bdf5ed7e43b760ac8b3b10c603b39309b580c3f538ee`.

Regenerate from the downloaded definitive run artifact, from this repository
root:

```sh
python tools/compact_rm3_moordyn_viz_corrected.py \
  --full-fine-dir /path/to/artifact/matlab-rm3-moordyn-viz-full-fine \
  --corrected-dir /path/to/artifact/matlab-rm3-moordyn-viz-corrected \
  --output tests/fixtures/matlab/rm3_moordyn_viz_corrected/source.npz
```

The hosted source workflow repeats the MATLAB run, regenerates this compact
archive, and compares every archived field, including scene signatures, with
the committed fixture. The portable native MoorDyn job runs the Python 80 s
motion and writes all 801 ParaView frames from its own body and line states.
The independent gate uses the unchanged published motion, PTO work, and
0.01 s fairlead limits plus full line path/load and scene signature limits.
For a local native run, set `WEC_SIM_MOORDYN_LIBRARY` to the pinned MoorDyn
shared library and run
`python -m pytest -q -s tests/test_rm3_moordyn_viz_corrected.py`.

## Measured limits

The source step refinement from 0.000625 to 0.0003125 s changed the largest
body position by 0.22515 mm, PTO force by 19.625 N, connection load by
170.76 N·m, 0.01 s fairlead tension by 60.99 N, line position by 0.683 mm,
and segment tension by 93.39 N. The corrected source's 0.1 s histories
matched the independently generated fine arm exactly. The original published
step differed from the corrected source by 0.20342 m in float surge and
104,189 N·m in connection pitch moment, so its strict expected failure is
still warranted.

The independent local Python run against the corrected scene had maximum
float surge/heave/pitch differences of 0.827 mm / 0.0733 mm / 7.89 µrad and
spar differences of 0.762 mm / 0.0772 mm / 7.89 µrad. PTO force differed by
177.1 N, and the largest signed PTO work error across the four published
windows was 0.00720%; connection pitch moment differed by 945.9 N·m.
Across 8,001 fairlead samples
the largest tension difference was 271.3 N. Across all three lines and 801
frames, the largest node and segment tension differences were 1.975 mm and
573.8 N. The Python-owned 801-frame scene matched the MATLAB wave sample and
extent signatures at five-decimal VTP precision; body centroid/bounds differed
by at most 0.84 mm, and line VTP positions by at most 1.98 mm.
