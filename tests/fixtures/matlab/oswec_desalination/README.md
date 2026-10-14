# Published OSWEC Desalination fixture

`source.npz` holds the seven CSV arrays used by the existing 300 s
Desalination network and coupled-motion gates. The arrays were loaded as
float64 from the successful [R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37853755507)
(`matlab-OSWEC_DESALINATION_SOURCE-baseline`, run commit
`5416034d09fbae9548a3f799f125381bf7352970`) and compressed without
rounding. The pinned WEC-Sim core is
`0753b2e47f2457c078751dcfe5d251d1767b80ab`; Applications is
`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`. The NPZ SHA-256 is
`211065d518990a9227a3e63ab17cfbbccff11dee4a0abc7173092ec1d047a7f9`.

| Source CSV suffix | SHA-256 |
| --- | --- |
| `simout1` | `0e8225c51a0754d7e874410c1f6ec47383a68fde00a6fa831acb5897d2638b47` |
| `relief` | `81de623dd4403a50da4b0e78dfc99a2014d80df8902afbb2e77e9e243976b2e8` |
| `cylinder` | `b6132435b1986fe319557cf92e75e6f4e57942d3417fbddbb0a7d36ae22c6a7d` |
| `valve1` | `7963a01276e353585c80d5bf2bf226c56b54625a6a30d53995068882ac494be2` |
| `components` | `406fd964571097a452a18654d5c07d4738e898b052e0de6eb78b080a5df47f9a` |
| `Desalination_body1` | `2ed16007b6b0afd67211f2c171b8ce534c6c14cac19acabbf2cefeba5c17b2cd` |
| `wave` | `1e78db924bb2b20cb9175d23a86141f2f9d20c5c83f1c62e7cbab894197a98fa` |

The coupled runner uses the existing 36-direction
[`oswec_multi_wave/oswec.h5`](../oswec_multi_wave/README.md) fixture. A local
R2025b regeneration from the pinned Applications `oswec.out` and `bemio.m`
matched all 65 comparable numerical HDF5 datasets; the largest absolute
difference was `1.46e-11` in a radiation kernel. The 10-direction
`examples/data/oswec.h5` changes excitation by as much as 2.28 N and cannot
serve as this published input.

The pinned Applications input sets `waves.phaseSeed = 1` for its 250-bin PM
sea. The portable coupled gate now generates its incident phases from that
published seed using WEC-Sim's MATLAB-compatible Threefry substream. All 250
phases are compared with the archived MATLAB values before the 300 s solve;
source phases, forces, and motion histories are comparison data, not inputs.

The tests retain their external CSV/HDF5 path for fresh MATLAB comparisons.
The portable fixture makes the already counted, full-duration Desalination
motion and hydraulic gates run in ordinary Python CI. It does not add a new
published-case badge claim. The source's raw chamber-pressure chatter remains
outside the paired force gate; the 0.1 s mean and net-work limits are unchanged.
