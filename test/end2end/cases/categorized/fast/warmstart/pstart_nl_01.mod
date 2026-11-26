# Primal start values without dual values (AMPL sends no dual values).
# Nonconvex: a local solver started at (2.1, -1.5) can end in another
# local minimum than from its default start (cold: -2.0839839).
# The suffix is declared here, so that it reads 0 when not reported.
suffix test_ini_pri OUT;

var x >= -3 <= 3 := 2.1;
var y >= -3 <= 3 := -1.5;

minimize f: sin(3*x) + cos(2*y) + 0.1*x*y;

s.t. c: x^2 + y^2 <= 8;
