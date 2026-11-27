# -------------------------------------------------------------
# inline_lin_const_01.mod
# (if c then 1) + (if not c then 1) is linearized to b + (1 - b):
# after inlining (cvt:pre:unnest bit 32) a constant without variables,
# inside exp; the second exp has a product term as well.
# Optimum -3 + e + exp(2) = 7.1073379274 at x = (0, 3, 0), brute force.
# -------------------------------------------------------------
var x{1..3} integer >= 0, <= 3;
minimize o: x[1] - x[2] + x[3]
  + exp((if x[1] >= 2 then 1) + (if not (x[1] >= 2) then 1))
  + exp(x[1]*x[2]/10 + (if x[2] >= 2 then 2) + (if not (x[2] >= 2) then 2));
s.t. c1: x[1] + x[2] + x[3] >= 2;
