# Exponential cone with a denominator of several terms:
# y*exp(z/(2*y+u)) <= x, the term u must not be dropped.
# Optimum y = 1, u = 3, exp(z/5) = 2.5, objective 3.4 - 2.5*ln(2.5).

var x;
var y >= 1 <= 5;
var z >= -5 <= 5;
var u >= 0 <= 3;
minimize o: x + 0.3*y - 0.5*z + 0.2*u;
s.t. c: y*exp(z/(2*y+u)) <= x;
