# -------------------------------------------------------------
# roundint_bnd_02.mod
# Rounded large integer bound used as big-M of an implication
# (cvt:pre:roundintbnd): x = 168000007 must stay feasible.
# AMPL presolve off: MP derives the bound.
# -------------------------------------------------------------
var x integer >= 0;
var b binary;
s.t. c: 0.1*x <= 16800000.7;
s.t. e: b == 0 ==> x <= 5;
maximize o: x;
