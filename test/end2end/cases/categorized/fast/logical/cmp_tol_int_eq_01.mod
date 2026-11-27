# -------------------------------------------------------------
# cmp_tol_int_eq_01.mod
# Integer x equal to p = 0.3/0.1 = 2.9999999999999996: within
# cvt:pre:inteps of 3, so x == p holds at x = 3.
# Optimum -3 at x = 3, y = 0.
# Exact (cvt:pre:inteps=-1 cvt:pre:ulps=0): x == p is false,
# y = 1, x = 5, optimum 5.
# -------------------------------------------------------------
param p := 0.3/0.1;

var x integer >= 0 <= 5;
var y binary;

minimize o: 10*y - x;

s.t. c: x == p or y = 1;
