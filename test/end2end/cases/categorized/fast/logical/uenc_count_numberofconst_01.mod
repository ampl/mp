set I default 1..3;
set DX default 1..9;

var x {I} in DX;
var y >=0 <=200;

## Check that even in positive context,
## we propagate mixed context (cvt:pre:ctx2count=0 default).
maximize Obj: y
              + count (x[1]==5, x[2]==3, x[3]==8);

s.t. LinCon: sum {i in I} x[i] == 15;

# Add more in negative context
s.t. LogCon: (x[1]==0 || x[3]==7) ==> y <= 10;

# This should add more negctx on x[3]==..
# and cause UEnc when ctx2count&2==0
s.t. CountCon: 2 <= numberof 4 in (x[1], x[2], x[3]);
