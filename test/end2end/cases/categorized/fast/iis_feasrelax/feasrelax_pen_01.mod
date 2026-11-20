# -------------------------------------------------------------
# feasrelax_pen_01.mod
# .rhspen semantics: c1 may not be relaxed (negative penalty),
# so c2 is relaxed: x + y = 2.
# -------------------------------------------------------------
suffix rhspen IN;
var x >= 0 <= 10;
var y >= 0 <= 10;
minimize o: x + y;
s.t. c1: x + y >= 2;
s.t. c2: x + y <= 1;
let c1.rhspen := -1;
let c2.rhspen := 1;
