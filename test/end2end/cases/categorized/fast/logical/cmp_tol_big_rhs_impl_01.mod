# -------------------------------------------------------------
# cmp_tol_big_rhs_impl_01.mod
# x > U never holds for x in [0, U] (U = 1e13), so y = 0 is
# not allowed.
# Optimum 1 at y = 1.
# -------------------------------------------------------------
param U default 1e13;

var x >= 0 <= U;
var y binary;

minimize o: y;

s.t. c: y = 0 ==> x > U;
