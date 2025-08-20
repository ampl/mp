set I default 1..3;
set DX default 1..9;

var x {I} in DX;
var y >=0 <=200;

## Propagation of exact/mixed context into inequalities
## (cvt:pre:ctx2ineq, cvt:pre:ctx2count).
maximize Obj: y
              + count (x[1]<=5, x[2]>=3, x[3]<=8)
              + x[1] - x[2] + x[3];

s.t. LinCon: sum {i in I} x[i] == 15;

# Add some opposite ones
s.t. LogCon: (x[1]<=5) ==> y <= 10;
