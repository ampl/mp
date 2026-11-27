# -------------------------------------------------------------
# round_bnd_fp_03.mod
# Nested rounding: the floor() bound as in round_bnd_fp_02
# inside ceil(min(...)/2). Optimum -11 at x = (1, -2, 5),
# where t = 2.
# -------------------------------------------------------------
var x1 integer >= 1 <= 1;
var x2 integer >= -2 <= 0;
var x3 integer >= 1 <= 5;

# ceil of half the smaller of two integer-valued terms
var t = ceil(min(max(2*x3, x1 - 2*x2),
                 floor((2*x3 - x2) / 3)) / 2);

minimize obj: -2*x1 + x2 - x3 - t;

s.t. c1: t >= 2;
