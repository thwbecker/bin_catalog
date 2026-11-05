# File formats

All files are whitespace separated ASCII, one record per line.
Longitudes and latitudes are in degrees, depths in km, angles in
degrees.

## Input catalogs

### AKI (default)

Ten columns per event:

```
lon lat depth strike dip rake mag plon plat tsec
```

`plon plat` are a second coordinate pair (for example projected or
longitude-shifted coordinates); they are read but then overwritten
with `lon lat` internally. `tsec` is the origin time in seconds since
1970-01-01 00:00 UTC. This last column is not part of the usual AKI &
Richards convention; build it with a UTC conversion, since a
local-time conversion offsets every event by the time zone. Events
with strike, dip, and rake all zero, or with a strike outside 0 to
360°, are reported and skipped. Moment is computed from the magnitude.
While reading, each event is compared with the previous 30 entries;
an event within 0.2 km horizontally, 0.4 km in depth, and 5 s in
origin time of an earlier one is treated as a duplicate and only the
larger magnitude is kept. The thresholds are `BC_DUP_HDIST_CLOSE`,
`BC_DUP_DEP_CLOSE`, and `BC_DUP_TIME_CLOSE` in `include/catalog.h`;
the number of duplicates found is reported on stderr.

### CMT (`-c`)

Thirteen columns per event, following the GMT `psmeca -Sm` layout
with the event title replaced by the UNIX time:

```
lon lat depth mrr mtt mff mrt mrf mtf exp plon plat tsec
```

Tensor components are in units of 10^exp dyn cm. The best fitting
double couple, scalar moment, and moment magnitude are derived on
input.

### ENG

Five columns, used by the magnitude utilities:

```
lon lat depth mag tsec
```

## Output of the summation programs

Six tensor components come first in every file, in the order

```
rr rt rp tt tp pp
```

in a spherical (r, θ, φ) system, with θ pointing south and φ east.
In the north, east, up system this corresponds to `zz -yz xz yy -xy xx`
(uu, -nu, eu, nn, -en, ee). The `psmeca -Sm` column order
`mrr mtt mff mrt mrf mtf` is obtained by picking columns 1, 4, 6, 2,
3, 5. The trace is removed from the summed tensors.

### `prefix.dx.dy.nmc.norm.dat`

Normalized Kostrov sum per bin or sample.

| Column | Content |
| --- | --- |
| 1-6 | normalized tensor (unit norm) |
| 7, 8 | lon, lat of the bin node or sample center |
| 9 | number of events |
| 10 | mean horizontal normal strain of the normalized tensor, (e1 + e2)/2 of the horizontal part |
| 11 | its standard deviation over the Monte Carlo realizations (nan without Monte Carlo) |

### `prefix.dx.dy.nmc.norm.sigma.dat` (only with `--monte-carlo` > 0)

Standard deviation of each of the six normalized components over the
realizations, followed by lon, lat, and the event count.

### `prefix.dx.dy.nmc.scaled.dat`

Moment-weighted sum, divided by the bin area and by the mean moment
per filled bin, times 10^6: columns 1-6 tensor, 7-8 lon lat, 9 event
count, 10 the corresponding mean horizontal normal strain.

### `prefix.dx.dy.nmc.clvd.dat`

`lon lat r_CLVD Gamma_CLVD` of the normalized tensor; both measures
are zero for a pure double couple.

## Output of the stress inversions

Rows are written only for bins with at least `--min-events` events
and a finite solution. Stress tensors are deviatoric and normalized;
the same component order as above applies.

### `prefix.dx.dy.s.dat` (Michael)

| Column | Content |
| --- | --- |
| 1-6 | stress tensor |
| 7, 8 | lon, lat |
| 9 | number of events |
| 10-15 | uncertainty of each component from random sweeps over the choice of nodal plane |
| 16 | norm of the uncertainty tensor |

### `prefix.dx.dy.smn.dat`

Difference between the normalized Michael stress and the normalized
Kostrov tensor (trace removed), same columns as `.s.dat`, with the
norm of the difference in the last column.

### `prefix.dx.dy.ds.dat` (Vavryčuk, fixed friction; `-F` ≥ 1)

Columns 1-6 stress, 7-8 lon lat, 9 event count, 10 mean fault
instability, 11 the friction used (0.6).

### `prefix.dx.dy.bs.dat` (Vavryčuk, best friction; `-F` ≥ 2)

Columns 1-6 stress, 7-8 lon lat, 9 event count, 10 mean instability
at the best friction, 11 best friction, 12 and 13 bootstrap mean and
standard deviation of the best friction (0 unless `-F 4`).

For `nsample_catalog` the file names are the same with `nmc` = 0; for
`psample_catalog` the summation files are `prefix.0.*.dat` and the
stress files `prefix.*.dat`.
