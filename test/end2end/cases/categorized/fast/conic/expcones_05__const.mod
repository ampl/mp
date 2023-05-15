# expcones_05__const.mod

var x {1..3};

s.t. Neg {i in 1..2}: x[i] <= 0;

maximize Obj:
   x[1] + x[2];

s.t. ExpCone01:
   -1.1 <= 3*x[2] * exp( 4*x[3] / (-3*x[2]) );

s.t. ExpCone02:
   -2*x[1] >= 0.3 * exp( 4*x[3] / (0.3) );

s.t. ExpCone03:
   1.1 >= -3*x[2] * exp( 0.4 / (-3*x[2]) );

s.t. LinCon:
   -x[1]-x[2]+x[3] == 1;
