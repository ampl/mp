## Check negated context propagation into arguments

var x {1..3} >= 0.0;

s.t. x3UB: x[3] <= 10;

s.t. ExpCone1: 1.1 >= x[3] -x[1] * exp(x[2]);

minimize Obj: x[1] + x[2] - x[3];  
