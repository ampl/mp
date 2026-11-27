# -------------------------------------------------------------
# cmp_bound_fp_01.mod
# A reified comparison must not be decided false from a computed
# bound that misses the right-hand side by a rounding error:
# the upper bound of (x + 1)/3 is computed as 1.9999999999999998,
# but (x + 1)/3 >= 2 holds at x = 5, where y >= 3 is then required.
# Optimum -8 at x = 4, y = 0 (x = 5 needs y >= 3: -7).
# -------------------------------------------------------------
var x integer >= 0 <= 5;
var y integer >= 0 <= 5;

minimize o: y - 2*x;

s.t. c: (x + 1) / 3 >= 2 ==> y >= 3;
