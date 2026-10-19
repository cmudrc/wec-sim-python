# Local R2025b WaveStar NMPC source

`source.npz` contains the 4,501 published 50 ms output samples from the full
225 s `WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive` run on Mali (macOS ARM64,
MATLAB R2025b). Pinned WEC-Sim is
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications is
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The temporary checkout
corrected only the application's `.stl` filename case to match Git's `.STL`
files. The source log reports WEC-Sim 7.1.0, 225 s end time, 0.05 s step, and
25 s wave ramp. `tests/matlab/reference_model_baseline.m` exports the source
body, PTO, wave, component, and controller signals.

The fixture selects time, float pitch and speed, PTO axial force and speed,
the 500 wave components, the saved wave, and the NMPC torque request from the
local CSV/MAT exports. SHA-256 of the compressed fixture is
`63b927e8d1208cac40406ed7feff750e65aeaeb6db987847104e7130b0ca184b`.
SHA-256 of the source `body.csv`, `pto.csv`, `components.csv`, `wave.csv`, and
`controller.mat` files, respectively:

- `8807d66659df5c0179441549c608d48a956921c20b3b029f13a52c6ac85c138d`
- `83f14a9dbd52bab6d4ac9134a528b55360417281c50a5d610f04bdc727e05e0e`
- `8e5b4b8e14501595f3e4da2dc424e32ab2629b9fe037177ca49e1ffb40376e37`
- `586ac5174bf559a8292cd3284a7ea4fda6787a6c1d67e54f0161ee6595bd4e2a`
- `cf62ed39db7ddb9f0c0e9730bb60a9989a796a2c04d013bf4c800f4117a4c4b0`

The regenerated local `wavestar.h5` is byte-identical to the pinned Python
fixture (`22702b94f052dfe6bb93d257006654d4e51246020a53f25b432d591bf8139d1a`).
The local and archived Linux MATLAB runs also have exactly identical 500 wave
components and 4,501 wave samples. Closed-loop pitch differs by up to
`0.00754261` rad and speed by `0.0502185` rad/s after NMPC activation;
requested torques first differ by more than `0.01` N m at 43 s. This fixture
supports a macOS ARM numerical comparison, not a platform-independent
physical parity claim. The case is qualified, while its strict badge remains
2/3 pending a portable full-duration physical gate.
