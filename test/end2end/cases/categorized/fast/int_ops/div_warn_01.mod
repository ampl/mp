# -------------------------------------------------------------
# div_warn_01.mod
# A divisor whose domain contains 0: warning divisor_domain_0
# (the divisor is treated as nonzero in the whole model).
# -------------------------------------------------------------
var x >= 0 <= 10;
var y integer >= -3 <= 3;
var b binary;
maximize o: if b then x/y else x;
s.t. c: x <= 4;
