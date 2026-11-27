# -------------------------------------------------------------
# feasrelax_ubpen_01.mod
# .ubpen of a variable: with x.ubpen < c.rhspen the bound x <= 1 is
# relaxed (x = 2); with a large x.ubpen, c is relaxed (x = 1).
# (CPLEX used penalties as FeasOpt preferences, i.e., inversely.)
# -------------------------------------------------------------
suffix ubpen IN;
suffix rhspen IN;
param pen_x default 0.5;
var x >= 0 <= 1;
var y >= 0 <= 10;
minimize o: y;
s.t. c: x + y >= 3;
s.t. d: y <= 1;
let x.ubpen := pen_x;
let c.rhspen := 1;
let d.rhspen := 10;
