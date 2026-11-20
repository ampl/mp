# -------------------------------------------------------------
# inline_quad_lin_01.mod
# Quadratic subexpression inside a linear expression defining another
# variable (cvt:pre:unnest bit 32), used in a nonlinear function
# -------------------------------------------------------------
var x >= 0, <= 2; var y >= 0, <= 2; var z integer >= 0, <= 3;
minimize o: -x - y + exp((x*y - x) + (if z >= 2 then 1) - (if not (z >= 2) then 1)) + z;
s.t. c1: x*y + z - (if x >= 1 then 1) <= 3;
