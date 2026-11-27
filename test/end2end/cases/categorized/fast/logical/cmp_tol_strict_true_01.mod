# -------------------------------------------------------------
# cmp_tol_strict_true_01.mod
# A strict comparison that holds exactly for every value of its
# expression is decided true: x > -5e-7 at x >= 0. Not left to the
# encoding with the gap cvt:mip:eps, which would require
# x >= 1e-4 - 5e-7 for y = 1.
# Optimum 0 at x = 0, y = 1.
# -------------------------------------------------------------
var x >= 0 <= 1;
var y binary;

minimize o: 1000*x + 10*(1 - y);

s.t. c: y = 1 ==> x > -5e-7;
