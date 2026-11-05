# Programs

All programs are built into `bin/` by `make`. Each prints a usage
message with `-h` or when called without arguments; the defaults
quoted below are those compiled in at the time of writing and are
printed by `-h`.

## Summation and stress inversion

The three main programs share most options and all output formats
(see [formats.md](formats.md)). They differ in how events are
assigned to output locations.

### `bin_catalog`

Regular binning. Events within ±dx/2 and ±dy/2 of a grid node are
assigned to that node; output coordinates are the node positions.

```
bin_catalog [options] catalog_file
```

| Option | Meaning | Default |
| --- | --- | --- |
| `-d, --dx val` | bin spacing in degrees (km with `-x`) | 0.2 |
| `-y, --dy val` | bin spacing in y, if different from dx | dx |
| `-m, --min-mag val`, `-M, --max-mag val` | magnitude range | 3.5, 6.5 |
| `-l, --min-lon`, `-r, --max-lon`, `-b, --min-lat`, `-t, --max-lat` | region | 232, 250, 30, 45 |
| `-Z, --min-depth val`, `-z, --max-depth val` | depth range, km | -10, 15 |
| `-n, --monte-carlo val` | number of Monte Carlo realizations; each rotates every mechanism by Gaussian random angles with σ = 30° about three axes and writes the spread to `.norm.sigma.dat` | 0 |
| `-w, --weighting val` | 0 none; 1 weight by temporal distance from the center of the catalog time span; 2 normalize by event frequency over time | 0 |
| `-p, --min-events val` | minimum number of events in a bin to attempt a stress inversion | 5 |
| `-F, --friction-solve val` | 0 Michael only; 1 add Vavryčuk at friction 0.6; 2 also scan friction 0 to 1; 3 scan 0.2 to 0.8; 4 scan 0 to 1 with bootstrap uncertainty | 2 |
| `-o, --out-prefix str` | output file prefix | kostrov |
| `-x, --xy` | treat coordinates as Cartesian x, y | off |
| `-c, --cmt` | read CMT format instead of AKI | off |

Output files are named `prefix.dx.dy.nmc.*.dat` for the summations
and `prefix.dx.dy.*.dat` for the stress tensors, with `nmc` the
number of Monte Carlo realizations.

### `nsample_catalog`

Same options as `bin_catalog`, but each grid node becomes a search
center. With `--nmin n` > 0, a node is kept if at least n events lie
within `--dist-max` km; with `--nmin -k`, the k nearest events are
used, provided they are all within `--dist-max`. Weighting `-w 1`
weights events by distance from the center. Defaults for `--nmin` and
`--dist-max` are 10 and 25 km. Output coordinates are the search
centers; file names follow `bin_catalog` with `nmc` = 0.

### `psample_catalog`

Like `nsample_catalog`, but the sample locations are read as
`lon lat` pairs from stdin and no grid is set up:

```
cat points.lonlat | psample_catalog [options] catalog_file
```

Output files are `prefix.0.*.dat` for the summations and
`prefix.*.dat` for the stress tensors; coordinates are the input
points.

### `solve_stress_one_bin`

Runs the Michael and Vavryčuk inversions on all mechanisms of one
catalog and prints the resulting tensors and friction statistics to
stderr; intended for checking a single bin. With `stdin_sdr` as the
file name, strike, dip, and rake triplets are read from stdin.

```
solve_stress_one_bin catalog.aki
echo "30 60 -90" | solve_stress_one_bin stdin_sdr
```

## Catalog utilities

| Program | Purpose |
| --- | --- |
| `merge_catalog a.aki b.aki out.aki` | merges two AKI catalogs and removes doubles, defined as entries within 1 km horizontally, 15 km in depth, and 5 s in origin time (`BC_MERGE_HDIST_CLOSE`, `BC_MERGE_DEP_CLOSE`, `BC_DUP_TIME_CLOSE` in `include/catalog.h`); of a double, the larger magnitude is kept |
| `m02dcfp catalog.cmt` | converts a CMT catalog to best fitting double couples, writing `lon lat depth strike dip rake strike2 dip2 rake2 mantissa exponent plon plat tsec` to stdout |
| `m02mag catalog.cmt` | writes `lon lat depth M0[Nm] Mw tsec` for each event of a CMT catalog |
| `calc_aux_plane` | reads `strike dip rake` from stdin and prints both nodal planes |
| `calc_gr [dM] [Mc] [mode]` | reads magnitudes from stdin and estimates the Gutenberg-Richter b value; modes 0 Aki maximum likelihood, 1 Utsu (1966)/Bender (1983) bin corrected with Shi & Bolt (1982) uncertainty (default), 2 Marzocchi, 3 b-positive, 4 all for comparison; defaults dM 0.01, Mc 2 |
| `calc_gr_time dt [dM] [Mc]` | reads `time magnitude` pairs (increasing time) from stdin and estimates b in windows of length dt, or in windows of N = -dt events when dt < 0 |

## Eigenvalue helpers

Small stand-alone drivers around the EISPACK routines, used by
`plot_kostrov` for Simpson style regime characterization:

| Program | Input on stdin | Output |
| --- | --- | --- |
| `eigen m n` | an m × n matrix | eigenvalues (Re, Im) and eigenvectors |
| `eigenvalues m n` | an m × n matrix | eigenvalues only |
| `eigen3ds` | rows of 3-D symmetric tensors as `11 12 13 22 23 33` | `e3 v3 e2 v2 e1 v1` with e1 > e2 > e3 |
| `eigenvalues3ds` | same | eigenvalues only |

## Scripts

| Script | Purpose |
| --- | --- |
| `scripts/plot_kostrov` | GMT 6 driver that runs `bin_catalog` or `nsample_catalog` for a region code from `region_parameters` and plots the summation or stress tensors as moment tensor symbols or strain crosses, optionally with topography, GPS velocities, plate boundaries, slab contours, and World Stress Map data. Positional arguments are documented at the top of the script; paths to external data default to `$HOME/data/` |
| `scripts/region_parameters` | prints region bounds, projection, annotation, and grid increments for numbered regions |
| `scripts/kostrov_to_lonlatazi1` | reduces a `.norm.dat` or `.s.dat` file to `lon lat azimuth` of the largest horizontal extensional eigenvector |
| `scripts/awk/*.awk` | helpers for tensor norms, horizontal principal axes, region strings, and World Stress Map sorting |

## Tests

`make test` (or `test/run_tests.sh [bin_dir]`) builds synthetic
catalogs and checks binning and sampling conventions, edge handling
(dateline, region boundaries, inverted ranges, invalid options), the
Vavryčuk inversion against reference values from the original MATLAB
code, reproducibility of the stress output, and termination of the
eigenvalue routines on non-finite input. `test/test_eigen.c` and
`test/test_weights.c` are unit tests built into `bin/`.
