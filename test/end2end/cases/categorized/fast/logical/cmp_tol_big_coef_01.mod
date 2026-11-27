# -------------------------------------------------------------
# cmp_tol_big_coef_01.mod
# Comparisons with large fractional coefficients (random model).
# In c2, -750000001.5*x2 == 1750000003.5 means x2 = -7/3, which
# an integer x2 cannot take: the bounds of the body miss the
# right-hand side by far, but 16 ulp of the right-hand side
# exceed cvt:pre:feastol. The converter decides it (before, it
# reached the solver, and Gurobi reported infeasible).
# c1: x2 = 2 or 3*x1 - 2*x2 < 11; c2: x1 > 0.
# Optimum -3 at x1 = 1, x2 = -2.
# -------------------------------------------------------------
var x1 integer >= -3 <= 4;
var x2 integer >= -2 <= 3;

minimize o: 3*x1 + 3*x2;

s.t. c1:
  -250000000.5*x2 == -500000001.0
  or 750000001.5*x1 - 500000001.0*x2 < 2750000005.5;

s.t. c2:                         # the premise always holds
  -750000001.5*x2 != 1750000003.5
  ==> -250000000.5*x1 < 0.0;
