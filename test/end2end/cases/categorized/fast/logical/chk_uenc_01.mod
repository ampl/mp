# -------------------------------------------------------------
# chk_uenc_01.mod
# Solution check bits 4 and 128 (intermediate expressions) with a
# unary encoding: no false warning for the encoding itself.
# Optimum -6 at x1 = -2, x2 = 0.
# -------------------------------------------------------------
var x1 integer >= -2 <= 5;
var x2 integer >= 0 <= 4;

minimize o: 3*x1 + 3*x2;

s.t. c1: (2*x2 - 3*x1 >= 6.5) or (-2*x2 != -4);
s.t. c2:
  atleast 1 (-2*x2 > -7.5,
             -3*x1 - 2*x2 != -11,
             -3*x1 >= 2.5);
s.t. c3: (2*x2 != 2) or (-x1 <= -2);
