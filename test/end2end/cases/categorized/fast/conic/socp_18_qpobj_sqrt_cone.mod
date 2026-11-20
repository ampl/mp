# Quadratic objective with coefficients 3 and 2 and a cone from sqrt(),
# which reaches the objective rewriting also on drivers that pass
# quadratic constraints natively. Optimum: x = 0.4, y = 0.6, objective 1.2.

var x >= 0 <= 10;
var y >= 0 <= 10;
minimize o: 3*x^2 + 2*y^2;
s.t. l: x + y >= 1;
s.t. c: sqrt(x^2 + y^2) <= 10;
