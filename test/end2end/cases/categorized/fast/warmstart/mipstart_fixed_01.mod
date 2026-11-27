# Integer variables all fixed by their bounds (with AMPL presolve off,
# the solver receives them): the model has no unfixed integer variables.
# Primal start values only.
var x integer >= 3 <= 3 := 3;
var y integer >= 2 <= 2 := 2;

minimize o: x + y;

s.t. c: x + y >= 1;
