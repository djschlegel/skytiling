#!/bin/bash
#SBATCH -J hsc_hir_anneal
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -c 256
#SBATCH -t 04:00:00
#SBATCH -o hsc_hir_anneal_%j.out
# Same as anneal_hsc_fibonacci.sh but starting from Atsushi's phase-1 layout
# (4 dithers on each of the 802 wide pointings), for a like-for-like comparison.
set -e
REPO=${SKYTILING_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}   # submit from the repo root
if [ ! -f "$REPO/bin/optimize_tiling_maps" ]; then echo "cannot find the repo at $REPO: submit from the repo root or set SKYTILING_REPO"; exit 1; fi
RUNDIR=${SKYTILING_RUNDIR:-$SCRATCH/skytiling/hsc_atsushi}
mkdir -p "$RUNDIR"; cd "$RUNDIR"
module load python
SKYTILING_PY=${SKYTILING_PY:-$REPO/../skytiling-venv/bin/python}
export NCPUS=${SLURM_CPUS_PER_TASK:-128}
NWORK=$(( NCPUS / 2 )); [ $NWORK -lt 1 ] && NWORK=1
DATA=$REPO/data/subaru
$SKYTILING_PY -u "$REPO/bin/optimize_tiling_maps" -t $DATA/hsc_mbq1_throughput.fits \
    --footprint $DATA/hsc_niji_wide_pointings.csv --footprint-radius 0.75 \
    --init-offsets "12.2,19.5 -19.5,12.2 -12.2,-19.5 19.5,-12.2" \
    --workers $NWORK --max-drift 0.75 --moves-per-tile 4 --delta-rot 10 --num-randoms 4000000 \
    --iters 300 --delta 0.1 --shrink 0.99 --seed 1 --plot --ra-center 180 --dec-center 0 --diameter 5 -o atsushi_
