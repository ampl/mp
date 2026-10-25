# Difference of PL terms in a defined variable,
# with the arguments fixed by constraints at and beyond the breakpoint.
# Needs AMPL presolve off and pl_linearize=0 to reach the solver.

var a;
var b;

s.t. FixA: a == 1;
s.t. FixB: b == 2;

var r = <<0; -1, 1>> a - <<0; -1, 1>> b;

maximize R: r;                                         # -1
