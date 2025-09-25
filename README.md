bin_catalog

Routines to read Aki, Engdahl, and Havard/Global CMT earthquake
catalogs and perform some basic computations, including computing
Gutenberg-Richter fits, binning, and merging catalogs.

Besides Kostrov summation (with some Monte Carlo error estimates), the
code can also compute Michael (1984) type best-fit stress tensors, as
well as Vavryčuk (2014) iterative inversions with fixed and best-fit
friction values.

The main Kostrov/stress inversion routines are bin_catalog for a
simple grid where each bin is summed and inverted for stress if there
are sufficient number of entries, as well as nsample_catalog which
determines the set of focal mechanisms to use based on a maximum
distance from point plus mininum number, or number of events within
distance criterion. See bin_catalog.c and nsample_catalog.c for the
options, plot_kostrov for some illustration as to how to use. 

Contains some versions of Andy Michael's code, as well as Vavryčuk
subroutines, converted to C, and GMT4.5.18 code - see their respective
copyright.

Michael, A.J., 1984. Determination of stress from slip data: Faults
and folds, J. Geophys. Res. 89, 11.517-11.526.

Michael, A.J., 1987. Use of focal mechanisms to determine stress: A
control study, J. Geophys. Res. 92, 357-368.

Vavryčuk, V., 2014. Iterative joint inversion for stress and fault
orientations from focal mechanisms, Geophysical Journal International,
199, 69-77.

(c) 2020 - 2025 by Thorsten Becker, thwbecker@post.harvard.edu


