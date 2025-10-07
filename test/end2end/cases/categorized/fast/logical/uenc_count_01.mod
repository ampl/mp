set I default 1..3;
set DX default 1..9;

var x {I} in DX;
var y >=0 <=200;

# Min!
minimize Obj: y
              + count (x[1]==5, x[2]==3, x[3]==8);

s.t. LinCon: sum {i in I} x[i] == 16;

# Add more in negative context
s.t. LogCon: (x[1]==6 || x[3]==17) ==> y >= 10;
