# Published RM3 hydraulic source fixture

`rm3.h5` and `source.npz` come from the successful [400 s R2025b MATLAB
run](https://github.com/cmudrc/wec-sim-python/actions/runs/37823137925)
(`matlab-RM3_HYDRAULIC_SOURCE-baseline`, artifact digest
`sha256:6d2d92fa081e8b4ea149eec72c27e9ab7096c5fbabcce037436f43b10d02c3b5`).
The workflow pins WEC-Sim core
`0753b2e47f2457c078751dcfe5d251d1767b80ab` and Applications
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The HDF5 is the exact
source file, SHA-256
`a1b357ef689bb12c7d958e23ab6036286865a729401264acb0a3c14a1fa40f44`.

The NPZ stores the source CSVs as float64 arrays at all 40,001 published
output times, plus all 500 incident-sea components. Only the two body arrays
were trimmed from 25 columns to their first 13: time, six positions, and six
velocities. All columns used by the existing source and trajectory gates are
retained. The NPZ SHA-256 is
`fd170a40c4857dcd78d1cbf267a00e5fe610673a3dd38bf3030d6b8764960d21`.
The original CSVs remain in the linked MATLAB artifact. Both test modules
still accept explicit fresh MATLAB output and HDF5 paths through their
environment variables.

These fixtures make the existing full-duration body, PTO, hydraulic, and
electrical gates run in ordinary Python CI with their numerical limits
unchanged. The case was already counted in the RM3 badge; this adds a durable
regression, not a new published-case gate.
