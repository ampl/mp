# -------------------------------------------------------------
# cmp_tol_bin_eq_01.mod
# 0.3*b == 0.1*3 for binary b: the right-hand side
# 0.30000000000000004 is within rounding error of the bound 0.3
# of 0.3*b, so the comparison holds at b = 1. Optimum 2 at
# b = z = 1. (AMPL evaluates the comparison exactly; z carries the
# converter's value.) Before, the reformulation had a big-M of
# 5.6e-17 and the solver's result depended on its tolerances.
# -------------------------------------------------------------
var b binary;
var z binary;

maximize o: b + z;

s.t. c: z = (if 0.3*b == 0.1*3 then 1 else 0);
