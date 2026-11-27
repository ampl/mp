# -------------------------------------------------------------
# cmp_tol_int_rhs_01.mod
# Integer x compared with p = 0.3/0.1 = 2.9999999999999996:
# within cvt:pre:inteps of 3, so x <= p allows x = 3, as it does
# as a bound or a row. Optimum 3 at x = 3, y = 0.
# Exact (cvt:pre:inteps=-1 cvt:pre:ulps=0): x <= 2, optimum 2.
# -------------------------------------------------------------
param p := 0.3/0.1;

var x integer >= 0 <= 5;
var y binary;

maximize o: x - 10*y;

s.t. c: x <= p or y = 1;
