# -------------------------------------------------------------
# cmp_tol_int_lt_01.mod
# Strict comparison of an integer variable with a fractional
# right-hand side: x < 3.5 allows x = 3, also when the right-hand
# side is not rounded first (cvt:pre:ineqrhs=0).
# Optimum 3.
# -------------------------------------------------------------
var x integer >= 0 <= 6;
var z binary;

maximize o: x;

s.t. c: x < 3.5 or z >= 1;
s.t. fz: z <= 0;
