# -------------------------------------------------------------
# cmp_tol_big_rhs_le_01.mod
# x <= U always holds for x in [0, U]. At U = 1e13 the
# floating-point margin of U (16 ulp = 0.03) exceeds
# cvt:pre:feastol; the converter decides with that margin
# instead of leaving b <==> x <= U to the solver, whose gap
# cvt:mip:eps is below the spacing of doubles at U.
# Optimum 1 at b = 1.
# -------------------------------------------------------------
param U default 1e13;

var x >= 0 <= U;
var b binary;

minimize o: b;

s.t. c: b <==> x <= U;
