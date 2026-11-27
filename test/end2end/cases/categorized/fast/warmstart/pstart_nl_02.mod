# Primal start values for a model without constraints (so there can be
# no dual values). From (2.1, -1.5), a local solver should end in the
# minimum -2.2512669 near (1.59, -1.61); from Gurobi's default start,
# in -2.0839839 near (-0.54, 1.58).
var x >= -3 <= 3 := 2.1;
var y >= -3 <= 3 := -1.5;

minimize f: sin(3*x) + cos(2*y) + 0.1*x*y;
