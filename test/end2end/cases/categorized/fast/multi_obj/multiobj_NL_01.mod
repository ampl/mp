#####################################################
# Emulated multi-objective pass whose first objective
# has a negative weight and a nonlinear term: the
# expression's relaxation direction must follow the
# direction the pass is solved in.
# Pass objective: -f1 + f2 = x - exp(x) + y, minimized:
# x = 2, y = 0 (f1 = exp(2), f2 = 2).
#####################################################

var x >= -2, <= 2;
var y >= 0, <= 3;

suffix objweight IN;

minimize f1: exp(x) + y suffix objweight -1;
minimize f2: x + 2*y suffix objweight 1;
