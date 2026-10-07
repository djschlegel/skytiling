#!/bin/bash
# Optimizer invocations used for the CFHT MegaCam solution.  Run from this directory.
# Stage 1: iterations 0-58, 10M randoms, 4 moves/tile (~13 min/iteration on a Mac)
#   Iter 58 | RMS = 0.6540 | sigma = 0.0564 deg
optimize_tiling \
      --force-fork \
      -t allsky-fib-165012-0.250.fits \
      -c cfht_40ccds.fits \
      --num-randoms 10000000 \
      --max-drift 1.50 \
      --delta 0.10 \
      --temp 0.0 \
      --iters 100 \
      --plot \
      --ra-center 10.0 \
      --dec-center 0.0 \
      --diameter 5.0 \
      --moves-per-tile 4 \
      --write-coverage

# Stage 2: resume from iteration 58 with 30M randoms (~25 min/iteration), through iteration 141
#   Initial RMS = 0.7977 (new randoms); Iter 59 | RMS = 0.7837 | sigma = 0.1000 deg
optimize_tiling \
      --force-fork \
      -t output_tiles_0058.fits \
      -c cfht_40ccds.fits \
      --num-randoms 30000000 \
      --max-drift 1.50 \
      --delta 0.10 \
      --temp 0.0 \
      --iters 100 \
      --plot \
      --ra-center 10.0 \
      --dec-center 0.0 \
      --diameter 5.0 \
      --moves-per-tile 4 \
      --write-coverage
