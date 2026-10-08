# Running on NERSC (Perlmutter)

`optimize_tiling_maps` parallelises over tiles the same way as
`optimize_tiling`: tiles are partitioned into layers of mutually distant
tiles (separation > 2 x (camera radius + max drift)), each layer's tiles are
distributed over forked worker processes that update a shared-memory
coverage array, and the layers are swept in turn.  One Perlmutter CPU node
(128 cores) is the right size; the HSC-Niji footprint run has ~3200 tiles
in ~50-100 layers, so 64 workers (one per physical core) are plenty.

## One-time setup

```
cd $HOME            # or wherever you keep code
git clone https://github.com/djschlegel/skytiling
module load python
python -m venv --system-site-packages $HOME/skytiling-venv
$HOME/skytiling-venv/bin/pip install numpy scipy astropy matplotlib
```

(The batch scripts create this venv next to the repo automatically if it is
missing; set `SKYTILING_PY` to use another python.)

## Submit

From the repo root (the scripts locate the repo through `SLURM_SUBMIT_DIR`;
to submit from elsewhere set `SKYTILING_REPO=/path/to/skytiling`):

```
sbatch nersc/anneal_hsc_fibonacci.sh        # Fibonacci start, 3208 tiles, 300 iterations
sbatch nersc/anneal_hsc_atsushi.sh          # Atsushi phase-1 start, same tile count
```

Output goes to `$SCRATCH/skytiling/hsc_fib/` (`fib_tiles_NNNN.fits`,
`fib_coverage_NNNN.png`, `fib_randoms_4000000.fits`) and the Slurm log
`hsc_fib_anneal_<jobid>.out` in the submission directory.  Edit `-A` for your
allocation and `-t` for the walltime; with 64 workers an iteration of the
3208-tile, 4M-random problem should take of order 5-10 s, so 300 iterations
is well under an hour.  Use `-q debug -t 00:30:00` for a quick test with
`--iters 20`.

## Resume

```
sbatch nersc/anneal_hsc_fibonacci.sh $SCRATCH/skytiling/hsc_fib/fib_tiles_0300.fits
```

continues from that checkpoint (the iteration counter continues from the
file name; `ra0, dec0` in the file keep the original drift reference) with a
smaller step.

## Tuning

* `--delta` / `--shrink`: proposal sigma [deg] and its per-iteration decay.
  The local runs used 0.1 x 0.96^n for 40 iterations and were still improving
  when sigma reached 0.02; the NERSC scripts use 0.1 x 0.99^n for 300.
* `--delta-rot`: sigma [deg] of the per-tile rotation-offset proposals (the
  four PAs stay 90 deg apart but the set turns as a whole); it shrinks with
  `--shrink` too.  0 freezes all tiles at `--rotations`.  `--init-rot random`
  starts from random offsets.
* `--moves-per-tile`: proposals per tile per iteration (best one is kept).
* `--max-drift`: cap on the distance from the starting position; it also sets
  the layer separation, so a larger value means fewer, bigger layers.
* `--temp`: Metropolis temperature (0 = greedy descent); energies are
  sum_f N Var(c_f) so a useful temperature is of order 1e-6 x that.
* `--use-good`: optimise the 0/1 footprint instead of the throughput maps.
