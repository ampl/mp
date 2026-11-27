# -------------------------------------------------------------
# ctx_ite_or_shared_cmp_01.mod
# The comparison x > 1 is used twice: as a disjunct of the root
# disjunction (positive context) and as the condition of an
# if-then-else whose result is bounded from above by a strict
# comparison (< 3). The condition needs both implications.
# Wrong before: x = 4 (violating the second constraint), -9.
# Correct optimum: x = 2, y = 0, objective -6.
# -------------------------------------------------------------

var x integer >= 0 <= 4;
var y integer >= 0 <= 4;

minimize o: -3*x + y;

s.t. c2: (x > 1) or (y > 2);
s.t. c3: (if x > 1 then x else 0) < 3;
