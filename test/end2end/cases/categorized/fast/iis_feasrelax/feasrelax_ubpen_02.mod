# -------------------------------------------------------------
# feasrelax_ubpen_02.mod
# .ubpen of a variable: with x.ubpen > c.rhspen, c is relaxed (x = 1);
# feasrelax_ubpen_01 relaxes the bound instead.
# (CPLEX used penalties as FeasOpt preferences, i.e., inversely.)
# -------------------------------------------------------------
suffix ubpen IN;
suffix rhspen IN;
param pen_x default 100;
var x >= 0 <= 1;
var y >= 0 <= 10;
minimize o: y;
s.t. c: x + y >= 3;
s.t. d: y <= 1;
let x.ubpen := pen_x;
let c.rhspen := 1;
let d.rhspen := 10;
