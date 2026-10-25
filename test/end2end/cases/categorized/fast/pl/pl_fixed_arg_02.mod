# PL terms in constraints whose (bounded) arguments are
# fixed by constraints, at and beyond the breakpoints.
# Needs AMPL presolve off and pl_linearize=0 to reach the solver.

set I = 1..4;
param a{i in I} =
    if 1==i then 2 else if 2==i then 0 else if 3==i then -1 else 0.5;

var x{I} >= -10, <= 10;
var y{I};

s.t. Fix{i in I}: x[i] == a[i];
s.t. AbsUB{i in I}: y[i] <= <<0; -1, 1>> x[i];        # nonconvex

maximize SumY: sum{i in I} y[i];                       # 3.5
