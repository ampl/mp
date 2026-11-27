# -------------------------------------------------------------
# cmp_tol_far_01.mod
# x <= 10, so x >= 1000 is false by a wide margin: the converter
# fixes b = 0 and passes no conditional comparison, also with
# a cvt:pre:feastol (1e-12) below the floating-point margin of
# 1000 (16 ulp = 1.8e-12).
# Optimum 10 at x = 10, b = 0.
# -------------------------------------------------------------
var x >= 0 <= 10;
var b binary;

maximize o: b + x;

s.t. c: b <==> x >= 1000;
