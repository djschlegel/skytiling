#!/bin/bash
#SBATCH -J hsc_fib_anneal
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -c 256
#SBATCH -t 04:00:00
#SBATCH -o hsc_fib_anneal_%j.out
#
# Anneal an HSC MBQ1 tiling of the HSC-Niji wide footprint on one Perlmutter
# CPU node, starting from a Fibonacci lattice with the same number of tiles
# as Hironao's phase-1 pattern (4 dithers x 802 pointings = 3208).
#
# Usage:  sbatch nersc/anneal_hsc_fibonacci.sh            (from the repo root)
#         sbatch nersc/anneal_hsc_fibonacci.sh fib_tiles_0100.fits   (resume)
# Set SKYTILING_PY to a python with numpy/scipy/astropy/matplotlib, or let
# the script create a venv next to the repo the first time (see nersc/README.md).

set -e
REPO=${SKYTILING_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}   # submit from the repo root
if [ ! -f "$REPO/bin/optimize_tiling_maps" ]; then echo "cannot find the repo at $REPO: submit from the repo root or set SKYTILING_REPO"; exit 1; fi
RUNDIR=${SKYTILING_RUNDIR:-$SCRATCH/skytiling/hsc_fib}
mkdir -p "$RUNDIR"; cd "$RUNDIR"
echo "repo $REPO, run dir $RUNDIR, $(date)"

module load python
if [ -z "$SKYTILING_PY" ]; then
    VENV=$REPO/../skytiling-venv
    if [ ! -x "$VENV/bin/python" ]; then
        echo "creating venv $VENV"
        python -m venv --system-site-packages "$VENV"
        "$VENV/bin/pip" install --quiet numpy scipy astropy matplotlib
    fi
    SKYTILING_PY=$VENV/bin/python
fi
export NCPUS=${SLURM_CPUS_PER_TASK:-128}
# physical cores only: the workers are numpy-bound, hyperthreads add little
NWORK=$(( NCPUS / 2 )); [ $NWORK -lt 1 ] && NWORK=1

DATA=$REPO/data/subaru
COMMON="-t $DATA/hsc_mbq1_throughput.fits --footprint $DATA/hsc_niji_wide_pointings.csv --footprint-radius 0.75
        --workers $NWORK --max-drift 0.75 --moves-per-tile 4 --delta-rot 10 --plot --ra-center 180 --dec-center 0 --diameter 5 -o fib_"

if [ -n "$1" ]; then
    # resume from a checkpoint (fib_tiles_NNNN.fits); reuse the randoms it was made with
    $SKYTILING_PY -u "$REPO/bin/optimize_tiling_maps" $COMMON --tiles "$1" --randoms fib_randoms_4000000.fits \
        --iters 200 --delta 0.03 --delta-rot 3 --shrink 0.99
else
    $SKYTILING_PY -u "$REPO/bin/optimize_tiling_maps" $COMMON --ntiles 3208 --margin 0.0 --num-randoms 4000000 \
        --iters 300 --delta 0.1 --shrink 0.99 --seed 1
fi
echo "done $(date)"
