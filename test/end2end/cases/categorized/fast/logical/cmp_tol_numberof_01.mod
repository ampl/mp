# -------------------------------------------------------------
# cmp_tol_numberof_01.mod
# numberof p in (...) with p = 0.3/0.1 = 2.9999999999999996
# counts integer variables at 3 (cvt:pre:inteps).
# Optimum 6: two variables at 3, one at 0.
# Exact (cvt:pre:inteps=-1 cvt:pre:ulps=0): no integer equals p,
# the model is infeasible.
# -------------------------------------------------------------
param p := 0.3/0.1;

var x {1..3} integer >= 0 <= 5;

minimize o: sum {i in 1..3} x[i];

s.t. c: numberof p in ({i in 1..3} x[i]) >= 2;
