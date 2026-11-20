# Example: g5sym_unrelaxed

As-built (unrelaxed) graphene | 5-AGNR | graphene junction: a 5-atom-wide
armchair ribbon of 10 units (H-passivated) lying AB-stacked at 3.40 A on two
graphene electrodes, with symmetric overlaps (7.88 / 8.59 A) and a 20 A gap
between the electrodes; 528 atoms (468 C, 60 H) in a 73.03 x 19.84 A cell
(`structure/`, a_CC = 1.4319 A). Transport along x, periodic along y.

```
./run_all.sh                    # the whole chain: kpath scf lo kcorr inputs qtbm gw
./run_all.sh inputs qtbm        # only some stages
MODEL=reference/lo/g5u_lo_kcorr_k3x3.npz ./run_all.sh inputs qtbm gw   # skip the DFT
```

`NP_SCF` (32), `NP_QTX` (64) and `NTHREADS` (32) set the MPI ranks of GPAW and
quatrex and the BLAS threads of the serial LO stages, `MPIFLAGS` extra
`mpirun` options (e.g. `--oversubscribe`, a hostfile) and `VENV` another
environment (default `../../.venv`). The number of quatrex ranks must divide
the 4096 energies.

| stage       | what                                                                   | time (reference run)   |
|-------------|------------------------------------------------------------------------|------------------------|
| `structure` | rebuild `structure/` (shipped; same options as the reference)          | seconds                |
| `kpath`     | band k path (Gamma-X, through K, Gamma-Y; 25 points each)              | seconds                |
| `scf`       | GPAW LCAO dzp PBE, 1 x 3 k (ky = 1/3 holds the folded Dirac point), 6384 orbitals | 25 min, 32 ranks |
| `lo`        | non-SCF 3 x 3 mesh, LCAO H(R), S(R) (5.9 GB npz), 468 p_z-like LOs     | ~1 h, 32 threads       |
| `kcorr`     | static correction on a 3 x 3 k mesh -> `g5u_lo_kcorr_k3x3.npz`          | minutes                |
| `inputs`    | quatrex device: 14 blocks (10 core + 2 x 2 padding), 820 LOs, E_F -1.024..+1.0235 eV, 4096 energies | minutes |
| `qtbm`      | ballistic T(E), ky = 3                                                 | 20 s, 64 ranks; 103 s, 8 ranks |
| `gw`        | NEGF GW, 50 SCBA iterations, mixing 0.5, eps_r 8                       | 1.6 h, 64 ranks        |

## Reference results (`reference/`)

- `lo/`: the LO model `g5u_lo_kcorr_k3x3.npz` (input of the `inputs` stage),
  `g5u_lo_summary.json`, `g5u_lo_kcorr_summary.json`, band and matrix plots.
  E_F = -4.0051 eV. Max band error of the LO model to the LCAO bands:

  | window      | static corr. at R = 0 | k3x3 correction |
  |-------------|-----------------------|-----------------|
  | +-0.25 eV   | 12.8 meV              | 3.2 meV         |
  | +-0.5 eV    | 23.9 meV              | 6.5 meV         |
  | +-1 eV      | 38.9 meV              | 16.4 meV        |

- `inputs/`: the quatrex inputs written from that model; check a new run with
  `../../.venv/bin/python ../../toolchain/check_inputs.py inputs reference/inputs`.
  Lead checks (`summary.json`): lead blocks periodic to 0.9 / 2.4 meV (left /
  right), couplings beyond 8 A drop at most 0.4 meV, band error of the cut 0.66 meV.
- `qtbm/transmission.png`: T(E) per ky and the ky average;
  `transmission.npz` the data (`energies` in eV absolute, `transmission_lr`
  per ky), from the release test run (8 ranks), which reproduced the plot.
- `gw/gw_ldos.png`, `gw/ldos_map.png`: SCBA convergence and ribbon / lead LDOS
  without and with GW; `quatrex_times.json` the profile. quatrex reports
  "SCBA did not converge after 49 iterations" at the iteration limit, while
  the maximum self-energy update has fallen to 3e-7.

The DFT stages write large files (`g5u_scf.gpw` 0.3 GB, `g5u_lo_lcao_HS.npz`
5.9 GB, `g5u_lo_lo_min_*.npz` 55 MB each) next to `run_all.sh`.
