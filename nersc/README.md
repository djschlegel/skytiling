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

## Re-anneal against a new throughput model

```
sbatch nersc/anneal_hsc_resume.sh data/subaru3 data/subaru2/anneal/nersc_atsushi2/tiles_0200.fits atsushi3_ 200
```

(those are the defaults, so `sbatch nersc/anneal_hsc_resume.sh` alone does the
same) copies the start file to `$SCRATCH/skytiling/hsc_atsushi3/atsushi3_tiles_0000.fits`
and anneals it against `data/subaru3/hsc_mbq1_throughput.fits` for 200
iterations with the smaller steps (`--delta 0.05 --delta-rot 5`), reusing the
4M randoms of the first run if `$SCRATCH/skytiling/hsc_fib/fib_randoms_4000000.fits`
exists.  The iteration counter restarts at 0 but `ra0, dec0` — and so the
0.75 deg drift cap from Atsushi's original positions — carry over.  About 20
minutes on one node.  When done, copy `atsushi3_tiles_0000.fits`,
`atsushi3_tiles_0200.fits` and the Slurm log's iteration lines into
`data/subaru3/anneal/nersc_atsushi3/`.

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

## Rectangular stripes (DESI Run 2 Niji stripes)

```
sbatch nersc/anneal_hsc_stripes.sh 1      # rms/mean objective
sbatch nersc/anneal_hsc_stripes.sh 2      # rms/mean^2: pulls tiles in harder
```

anneals a lattice + dither start (2548 tiles) over `data/subaru3/niji_stripes_2p8_hsc.csv`
(the Niji stripes with their RA bounds trimmed to the HSC-Niji wide footprint; a fourth
argument names another rectangle file, e.g. `niji_stripes_2p8.csv` for the untrimmed DESI bounds)
against the `data/subaru3` model, 300 iterations (`--objective cv`, see
`data/subaru3/README.md`); output `$SCRATCH/skytiling/hsc_stripes/stripes_p<p>_tiles_NNNN.fits`
with a `pid` (pointing id) column.  Both runs share the randoms file written by the first.
