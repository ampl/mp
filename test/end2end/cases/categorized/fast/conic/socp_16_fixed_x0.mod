# A cone whose first variable is fixed by its bounds (not the shared
# constant 1): sqrt(x^2+y^2) <= 3 as a quadratic constraint must be
# x^2+y^2 <= 9, not <= 3. Optimum x = y = 3/sqrt(2), objective 3*sqrt(2).
# Needs presolve 0, else AMPL substitutes t.

var x >= -10 <= 10;
var y >= -10 <= 10;
var t >= 3 <= 3;
maximize o: x + y;
s.t. s: sqrt(x^2 + y^2) <= t;
