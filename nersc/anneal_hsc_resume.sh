#!/bin/bash
#SBATCH -J hsc_resume
#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -c 256
#SBATCH -t 02:00:00
#SBATCH -o hsc_resume_%j.out
# Re-anneal an existing HSC tiling against a (new) throughput model, keeping
# the original drift reference (ra0, dec0) of the starting file.  This is how
# the data/subaru2 solution was made from the data/subaru one, and the
# data/subaru3 solution from the data/subaru2 one.
#
#   sbatch nersc/anneal_hsc_resume.sh [MODEL_DIR] [START_TILES] [PREFIX] [ITERS]
#
# defaults: MODEL_DIR=data/subaru3, START_TILES=data/subaru2/anneal/nersc_atsushi2/tiles_0200.fits,
#           PREFIX=atsushi3_, ITERS=200; run directory $SCRATCH/skytiling/hsc_<PREFIX without trailing _>
#
# The start file is copied to <PREFIX>tiles_0000.fits in the run directory, so
# the iteration counter restarts at 0 and the output is <PREFIX>tiles_0000.fits
# (= the start), ..., <PREFIX>tiles_<ITERS>.fits.  To continue an interrupted
# run, resubmit with START_TILES = the last checkpoint written (its counter is
# then continued instead).
set -e
REPO=${SKYTILING_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}   # submit from the repo root
if [ ! -f "$REPO/bin/optimize_tiling_maps" ]; then echo "cannot find the repo at $REPO: submit from the repo root or set SKYTILING_REPO"; exit 1; fi
MODEL=${1:-$REPO/data/subaru3}; [ -d "$MODEL" ] || MODEL=$REPO/$MODEL
START=${2:-$REPO/data/subaru2/anneal/nersc_atsushi2/tiles_0200.fits}; [ -f "$START" ] || START=$REPO/$START
PREFIX=${3:-atsushi3_}
ITERS=${4:-200}
RUNDIR=${SKYTILING_RUNDIR:-$SCRATCH/skytiling/hsc_${PREFIX%_}}
mkdir -p "$RUNDIR"
module load python
SKYTILING_PY=${SKYTILING_PY:-$(which python)}
export NCPUS=${SLURM_CPUS_PER_TASK:-128}
NWORK=$(( NCPUS / 2 )); [ $NWORK -lt 1 ] && NWORK=1
# reuse the 4M randoms of the first HSC run if they exist, otherwise generate them here
RANDOMS=$SCRATCH/skytiling/hsc_fib/fib_randoms_4000000.fits
case "$START" in
  "$RUNDIR"/${PREFIX}tiles_[0-9][0-9][0-9][0-9].fits) TILES=$START ;;            # continuing this run
  *) TILES=$RUNDIR/${PREFIX}tiles_0000.fits; cp "$START" "$TILES" ;;
esac
cd "$RUNDIR"
if [ -f "$RANDOMS" ]; then RANDARG="--randoms $RANDOMS"; else RANDARG="--num-randoms 4000000"; fi
$SKYTILING_PY -u "$REPO/bin/optimize_tiling_maps" -t "$MODEL/hsc_mbq1_throughput.fits" \
    --footprint "$REPO/data/subaru/hsc_niji_wide_pointings.csv" --footprint-radius 0.75 $RANDARG \
    --tiles "$TILES" --workers $NWORK --max-drift 0.75 --moves-per-tile 4 \
    --delta 0.05 --delta-rot 5 --shrink 0.99 --seed 3 --iters $ITERS \
    --plot --ra-center 180 --dec-center 0 --diameter 5 -o "$PREFIX"
