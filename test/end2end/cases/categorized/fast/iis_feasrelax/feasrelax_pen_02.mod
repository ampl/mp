# -------------------------------------------------------------
# feasrelax_pen_02.mod
# .rhspen semantics: relaxing c1 is free (zero penalty),
# so c1 is relaxed: x + y = 0.
# -------------------------------------------------------------
suffix rhspen IN;
var x >= 0 <= 10;
var y >= 0 <= 10;
minimize o: x + y;
s.t. c1: x + y >= 2;
s.t. c2: x + y <= 1;
let c1.rhspen := 0;
let c2.rhspen := 1;
