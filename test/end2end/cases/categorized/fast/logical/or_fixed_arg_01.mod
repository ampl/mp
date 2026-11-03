# OR/AND constraints with an argument fixed to 1 by its bounds.
# With AMPL and solver presolve off, the fixed argument reaches
# native AND/OR constraints (Xpress requires binary arguments).
var x {1..3} binary;
var c binary >= 1;

minimize Obj: x[1] + 2*x[2] + 3*x[3];

s.t. Or1: x[1] = 1 or x[2] = 1 or c = 1;
s.t. Or2: (x[2] = 1 and c = 1) or x[3] = 1;
