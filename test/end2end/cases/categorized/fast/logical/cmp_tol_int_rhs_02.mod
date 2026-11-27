# -------------------------------------------------------------
# cmp_tol_int_rhs_02.mod
# Integer x, x <= 3 - 1e-7 in a disjunction: the right-hand side
# within cvt:pre:inteps of 3 is 3, as for the bound or the row
# (where the solver accepts x = 3 by its own tolerance).
# Optimum 3. Exact (cvt:pre:inteps=-1 cvt:pre:ulps=0): 2.
# -------------------------------------------------------------
var x integer >= 0 <= 6;
var z binary;

maximize o: x;

s.t. c: x <= 3 - 1e-7 or z >= 1;
s.t. fz: z <= 0;
