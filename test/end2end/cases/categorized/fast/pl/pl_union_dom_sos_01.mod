## A nonconvex PL of a variable with an interval-union domain.
## With pl_linearize=1, AMPL's SOS2 binaries for the PL term precede the
## SOS1 variables for the domain; ASL drivers that remove the SOS2
## binaries must renumber the members of both sets.
## The PL term is in a constraint, not the objective, so that AMPL
## versions before 20261027 also give the right model.

var x in {1} union interval[2.5,4] union {6} union interval[8,9];
var t;
s.t. C: t >= <<1.5,3,6; -2,2,-2,2>> x;
minimize O: t;
