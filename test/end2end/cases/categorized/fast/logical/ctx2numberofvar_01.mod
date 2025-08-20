set I default 1..3;
set DX default 1..9;

var x {I} in DX;
var y >=0 <=200;

## Propagation of exact/mixed context into inequalities
## (cvt:pre:ctx2ineq, cvt:pre:ctx2count).
maximize Obj: y
              + x[1] - x[2] + x[3];

s.t. LinCon: sum {i in I} x[i] == 15;

s.t. NumberofVar: 2
                  <= numberof y+10*x[1]
                              in (x[1] + x[2]^2, x[2] + x[3]^3, x[3]^2-x[1]^2);
