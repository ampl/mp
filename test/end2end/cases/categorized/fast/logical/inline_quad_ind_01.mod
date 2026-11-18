# -------------------------------------------------------------
# inline_quad_ind_01.mod
# Quadratic subexpression (a product) in an indicator body, with a
# linearized NOT next to it: inlining makes a quadratic indicator
# -------------------------------------------------------------
var x >= 0, <= 3; var y >= 0, <= 3; var z integer >= 0, <= 3; var b binary;
minimize o: -x - y - z + 3*b;
s.t. c1: b ==> x*y + (if not (z <= 1) then 1) <= 2;
s.t. c2: (not b) ==> x + y + z <= 2;
