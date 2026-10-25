# PL terms in the objective whose (unbounded) arguments are
# fixed by constraints, at and beyond the breakpoints.
# Needs AMPL presolve off and pl_linearize=0 to reach the solver.
# Fixing at a breakpoint used to break the SOS2 -> MIP conversion.

set I = 1..4;
param a{i in I} =
    if 1==i then 2 else if 2==i then 0 else if 3==i then -1 else 0.5;

var x{I};

s.t. Fix{i in I}: x[i] == a[i];

maximize SumAbs: sum{i in I} <<0; -1, 1>> x[i];       # 3.5
