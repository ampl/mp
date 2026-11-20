# -------------------------------------------------------------
# int_ops_03.mod
# x div y with a variable divisor whose abs(y) appears elsewhere
# (in the objective, positive context only): the divisor's bound
# |rem| <= |y| - eps reused that abs and stopped with "Part of
# expression redefinition could be lost". Optimum 10 at y = -1.
# -------------------------------------------------------------
var x integer >= 0 <= 9;
var y integer >= -2 <= 2;
var q integer >= -20 <= 20;
minimize o: 10*abs(y) + q + x;
s.t. c: x >= 3;
s.t. g: q = x div y;
