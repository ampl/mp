# -------------------------------------------------------------
# cmp_tol_eq_scaled_01.mod
# 1e-3*x == 1e-3*(3 + 1e-7), x integer: x == 3.0000001 after
# division by the coefficient, within cvt:pre:inteps of 3.
# The disjunction is true at the root: x is fixed at 3 (before:
# contradictory bounds [4, 3], infeasible).
# Optimum 3.
# -------------------------------------------------------------
var x integer >= 0 <= 6;
var z binary;

maximize o: x;

s.t. c: 1e-3*x == 1e-3*(3 + 1e-7) or z >= 1;
s.t. fz: z <= 0;
