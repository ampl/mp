# -------------------------------------------------------------
# inline_quad_in_lin_01.mod
# A quadratic subexpression (x*y) inside the linear expression defining
# the if-then-else result, within a nonlinear function.
# -------------------------------------------------------------
var x >= 0, <= 2;
var y >= 0, <= 2;
var b binary;
minimize o: -x - y - b + sqrt(1 + (if b == 1 then x*y + x - 1 else x - y)^2);
s.t. c1: x + y <= 3;
