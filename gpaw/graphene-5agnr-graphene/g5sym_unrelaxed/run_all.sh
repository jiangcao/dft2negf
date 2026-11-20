#!/usr/bin/env bash
# g5sym_unrelaxed: GPAW -> LO model -> quatrex device -> quatrex (QTBM, GW) for the
# as-built (unrelaxed) symmetric AB-stacked graphene | 5-AGNR (10 units) | graphene
# junction, 528 atoms (468 C, 60 H) in a 73.0 x 19.8 A cell (structure/).
#
# Stages (run in this order; give some of them as arguments to run only those):
#   structure  rebuild structure/ with toolchain/build_gnr_junction.py (shipped already)
#   kpath      k path of the band lines (bands_junction.py --kpath-only: the reference
#              bands are then taken from the LCAO H(R) itself, lo_junction.py --ref-tb)
#   scf        GPAW LCAO dzp PBE SCF on a 1 x 3 k mesh, NP_SCF MPI ranks (default 32)
#   lo         LCAO H(R), S(R) on a 3 x 3 mesh, minimal LO model (1 p_z-like LO per C);
#              serial, NTHREADS BLAS threads (default 32); writes a 5.9 GB g5u_lo_lcao_HS.npz
#   kcorr      static hybridization correction of the LO model on a 3 x 3 k mesh
#              -> g5u_lo_kcorr_k3x3.npz, and pictures of the LO matrices
#   inputs     quatrex device: inputs/{hamiltonian,overlap,coulomb_matrix}.h5, ...,
#              qtbm/quatrex_config.toml, gw/quatrex_config.toml
#   qtbm       quatrex ballistic transmission (NP_QTX ranks, default 64), qtbm/transmission.png
#   gw         quatrex NEGF GW (SCBA, 50 iterations, ~1.6 h on 64 ranks), gw/*.png
#
#   ./run_all.sh                         # kpath .. gw (the full chain, a few hours)
#   ./run_all.sh inputs qtbm             # device and QTBM only, from the shipped
#   MODEL=reference/lo/g5u_lo_kcorr_k3x3.npz ./run_all.sh inputs qtbm    # LO model
#
# The DFT stages skip themselves when their output exists (delete it to redo).
# Results of the release run are in reference/; ../../toolchain/check_inputs.py
# compares a new inputs/ with reference/inputs/.
set -euo pipefail
cd "$(dirname "$0")"
EX=$(pwd)
ROOT=$(cd ../.. && pwd)
TC=$ROOT/toolchain
VENV=${VENV:-$ROOT/.venv}
PY=$VENV/bin/python
NP_SCF=${NP_SCF:-32}
NP_QTX=${NP_QTX:-64}
NTHREADS=${NTHREADS:-32}
# extra mpirun options, e.g. MPIFLAGS="--oversubscribe" or a hostfile
MODEL=${MODEL:-g5u_lo_kcorr_k3x3.npz}
S=structure/graphene_5AGNR_graphene_u10_sym.traj

threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
mpi() { mpirun -np "$1" ${MPIFLAGS:-} -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS -x MKL_NUM_THREADS "${@:2}"; }
say() { echo "[$(date '+%F %T')] $*"; }

stage_structure() {
    # all builder options spelled out: a_CC 1.4319 A (LCAO dzp PBE + D3(BJ) graphene),
    # interlayer distance 3.40 A, symmetric AB overlaps, 4 lead periods beyond each ribbon end
    (cd structure && "$PY" "$TC/build_gnr_junction.py" \
        --width 5 --units 10 --gap 20 --lead 17 --ny 8 --d 3.40 --vacuum 15 \
        --acc 1.4319 --symmetric --prefix graphene_5AGNR_graphene_u10_sym)
}

stage_kpath() {
    [ -f g5u_bands.npz ] || "$PY" "$TC/bands_junction.py" $S --prefix g5u_bands \
        --kpts 1 3 1 --npts 25 --kpath-only
}

stage_scf() {
    # --kpt-par 1: all ranks per k (complex pzheevd crashed on 8-rank groups)
    [ -f g5u_scf.gpw ] && return
    threads 1
    mpi "$NP_SCF" "$VENV/bin/gpaw" python "$TC/bands_junction.py" -- \
        $S --prefix g5u --kpts 1 3 1 --kpt-par 1 --scf-only > g5u_scf.out 2>&1
}

stage_lo() {
    threads "$NTHREADS"
    [ -f g5u_lo_lo_min_zEF.npz ] || "$PY" "$TC/lo_junction.py" g5u_bands --gpw g5u_scf.gpw \
        --ref-tb --no-ext --prefix g5u_lo --emax 1.0 > g5u_lo.out 2>&1
}

stage_kcorr() {
    threads "$NTHREADS"
    "$PY" "$TC/lo_kmesh_correction.py" g5u_lo_refbands g5u_lo --prefix g5u_lo_kcorr \
        --meshes 3 3 > g5u_lo_kcorr.out 2>&1
    "$PY" "$TC/plot_lo_matrices.py" g5u_lo_kcorr_k3x3.npz $S --prefix g5u_lo_k3x3 >> g5u_lo_kcorr.out 2>&1
}

stage_inputs() {
    # core of 10 blocks of 2 lead periods; QTBM / GW grid E_F -1.024 .. +1.0235 eV, 4096
    # points (0.5 meV); GW: 50 SCBA iterations, mixing 0.5, eps_r 8
    threads 8
    "$PY" "$TC/make_quatrex_inputs.py" "$MODEL" $S --out . \
        --gw-iter 50 --mixing 0.5 --epsr 8 --ewin -1.024 1.0235 --ne 4096 > make_inputs.out 2>&1
}

quatrex() {   # $1 = qtbm or gw
    local ne
    ne=$(sed -n "s/^energy_window_num = \([0-9]*\).*/\1/p" "$1/quatrex_config.toml")
    (( ne % NP_QTX == 0 )) || { echo "energy_window_num $ne not divisible by $NP_QTX ranks" >&2; exit 1; }
    threads 1
    (cd "$1" && mpi "$NP_QTX" "$VENV/bin/quatrex" run quatrex_config.toml > "$1.out" 2>&1)
}

stage_qtbm() { quatrex qtbm; (cd qtbm && "$PY" plot_T.py); }
stage_gw() { quatrex gw; (cd gw && "$PY" plot_gw.py && "$PY" plot_ldos_map.py); }

ALL=(kpath scf lo kcorr inputs qtbm gw)
STAGES=("$@")
[ ${#STAGES[@]} -gt 0 ] || STAGES=("${ALL[@]}")
for st in "${STAGES[@]}"; do
    declare -F "stage_$st" > /dev/null || { echo "unknown stage '$st' (structure ${ALL[*]})" >&2; exit 1; }
done
for st in "${STAGES[@]}"; do
    say "stage $st"
    "stage_$st"
done
say "done: ${STAGES[*]}"
