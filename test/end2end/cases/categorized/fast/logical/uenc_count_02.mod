set I default 1..3;
set DX default 1..9;

var x {I} in DX;
var y >=0 <=200;

## Check that even in positive context,
## we propagate mixed context (cvt:pre:ctx2count=1 default).
maximize Obj: y
              + count (x[1]==5, x[2]==3, x[3]==8);

s.t. LinCon: sum {i in I} x[i] == 15;

# Add more in negative context
s.t. LogCon: (x[1]==0 || x[3]==7) ==> y <= 10;
