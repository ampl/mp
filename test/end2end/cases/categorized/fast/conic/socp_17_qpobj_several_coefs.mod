# Several squared objective terms with different coefficients, linear terms
# and a constant, plus a quadratic constraint recognized as a cone: a
# mismatch between coefficients and variables or a wrong scaling shows up.
# Optimum: x = 1, y = 2, objective 3.

var x >= 0 <= 10;
var y >= 0 <= 10;
minimize o: 2*x^2 + 0.5*y^2 - 4*x - 2*y + 7;
s.t. c: x^2 + y^2 <= 100;
