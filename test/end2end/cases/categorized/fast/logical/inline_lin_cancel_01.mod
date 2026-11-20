# -------------------------------------------------------------
# inline_lin_cancel_01.mod
# Nested linearizations of NOT whose terms cancel when inlined
# (cvt:pre:unnest bit 32), e.g. not(not(.)): a linear expression
# can become constant. Optimum -2 by brute force.
# -------------------------------------------------------------
var x{1..2} integer >= 0, <= 3;
var b binary;
minimize o: b + x[1] - x[2];
s.t. c1: ((x[1] <= 1) <==> not (x[2] >= 1)) or (b >= 1);
s.t. c2: (x[1] + x[2] >= 2) <==> not (not (x[1] + x[2] >= 2));
