# A cone whose first variable is fixed by its bounds at a negative value:
# sqrt(x^2+y^2) <= t with t = -3 is infeasible. As a quadratic constraint
# the right-hand side must keep the sign (x^2+y^2 <= -9), not become 9.
# Needs presolve 0, else AMPL substitutes t.

var x >= -10 <= 10;
var y >= -10 <= 10;
var t >= -3 <= -3;
maximize o: x + y;
s.t. s: sqrt(x^2 + y^2) <= t;
