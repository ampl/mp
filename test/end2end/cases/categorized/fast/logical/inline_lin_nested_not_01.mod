# -------------------------------------------------------------
# inline_lin_nested_not_01.mod
# Linearized NOT inside equivalences: the linear functional constraints
# of the equivalences contain those of the NOTs (cvt:pre:unnest bit 32).
# -------------------------------------------------------------
var x{1..3} integer >= 0, <= 3;
minimize o: x[1] - 2*x[2] + x[3];
s.t. c1: (x[1] + x[2] >= 3) <==> not (x[3] <= 1);
s.t. c2: not (x[2] >= 2) <==> not (not (x[1] <= 1));
