#############################################
## Test inlining of nonlinear expressions
## in constraint and objectives
#############################################

var x {i in 1..4} in interval[-3, 4] := 3 * (-1)^i;

var sin1 = sin(x[1] + 3*x[3]);
var log1 = log(x[3]^2);
var log2 = log( tan(2*x[2] + 4*x[4] - 8) );

minimize Obj:
   x[1] - x[2] + x[3] - x[4]
   + log1
   + log2
   - sin1;

s.t. C1:
   -sin1
   + log2
   + cos(x[3] - 12*x[4]^3)
   >= 5;

s.t. Fix_Log1: log1 == 1.5;

s.t. Bound_Sin1: sin1 >= 0.5;

s.t. Bound_Log2: log2 <= 6.2;

let x[1] := 0.5;
let x[2] := 4;
let x[3] := -2;
let x[4] := 3.5;

display _obj;
