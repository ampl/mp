# expcones_04__negvar.mod

var x {1..3};

s.t. Neg {i in 1..2}: x[i] <= 0;

maximize Obj:
   x[1] + x[2];

s.t. ExpCone:
   -2*x[1] >= -3*x[2] * exp( 4*x[3] / (-3*x[2]) );

s.t. LinCon:
   -x[1]-x[2]+x[3] == 1;
