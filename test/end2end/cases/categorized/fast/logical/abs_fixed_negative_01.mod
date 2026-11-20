# -------------------------------------------------------------
# abs_fixed_negative_01.mod
# abs() of a variable fixed at a negative value is a constant.
# Wrong before: the constant was taken for a variable index, the
# comparison abs(x) > 1 was decided false, "infeasible".
# Correct: x = -2, objective -2.
# -------------------------------------------------------------

var x integer >= -2 <= -2;

minimize o: x;

s.t. c1: abs(x) > 1;
