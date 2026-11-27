# -------------------------------------------------------------
# cmp_tol_neg_scaled_01.mod
# not (1e-5*x > 5e-5) is 1e-5*x <= 5e-5, which holds at x = 5;
# it must not become 1e-5*x == 0 near the bound 0.
# Optimum 5.
# -------------------------------------------------------------
var x integer >= 0 <= 8;

maximize o: x;

s.t. c: not (1e-5*x > 5e-5);
