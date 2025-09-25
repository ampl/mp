###########################################
## Test unification of inequalities
###########################################

var b {1..3} binary;

var x {1..5} integer >=0 <=5;

minimize Obj: b[1] + b[2] - x[2] - 2*x[4] + x[1] - 5*x[3] - x[5];

s.t. C1: b[1]>0.5 ==> x[2] - x[4] >= -2;
s.t. C2: b[2]<0.5 ==> x[2] - x[4] < -2;
s.t. C3: x[2] - x[4] > -3 || x[5] < 4;
s.t. C4: x[2] - x[4] <= -3 ==> x[5] <= 3;

## Unify (ineq)-LB to (==LB)
s.t. C5: exactly 1 (x[1] - 3*x[3] > -15, x[5] > 4, x[5] < 1);

## See if this reuses the above information that it is "exactly 1"
s.t. C6: exactly 1 (x[1] - 3*x[3] >= -14, x[5] == 5, x[5] == 0)
           ==> x[3] > 2;
