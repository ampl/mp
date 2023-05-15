# expcones_01__plain.mod
# From https://docs.mosek.com/latest/capi/tutorial-ceo-shared.html

var x {1..3};

s.t. Nonneg {i in 1..2}: x[i] >= 0;

minimize Obj:
   x[1] + x[2];

s.t. ExpCone:
   x[1] >= x[2] * exp( x[3] / x[2] );

s.t. LinCon:
   sum {i in 1..3} x[i] == 1;
