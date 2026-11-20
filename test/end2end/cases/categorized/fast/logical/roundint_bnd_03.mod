# -------------------------------------------------------------
# roundint_bnd_03.mod
# Infeasible: dl gives x >= 599999999.5, i.e. x >= 600000000 for
# integer x; c0 then requires x < 600000000. A fractional data bound
# must not be rounded to the nearest integer (cvt:pre:roundintbnd).
# AMPL presolve off: MP derives the bound.
# -------------------------------------------------------------
var x integer >= 599999950 <= 600000050;
s.t. dl: 3*x >= 1799999998.5;
s.t. c0: (x >= 600000000) ==> (x < 600000000);
minimize o: x;
