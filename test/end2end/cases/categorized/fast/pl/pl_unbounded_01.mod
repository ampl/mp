# Nonconvex PL term whose argument has no bounds:
# x in [-4, 2] is implied only by 2-variable constraints.
# Needs AMPL presolve off and pl_linearize=0 to reach the solver.
# The MIP reformulation of SOS2 needs a domain for x:
# without cvt:pl:domain, an unbounded argument is an error.

var x;
var w >= 0 <= 1;

s.t. c1: x + w <= 2;
s.t. c2: x - w >= -4;

maximize AbsX: <<0; -1, 1>> x;                         # 4
