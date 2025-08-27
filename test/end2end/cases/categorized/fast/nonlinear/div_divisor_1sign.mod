## Division where divisor is 1-sign.
## Check we don't introduce indicators with acc:div=0.

var x >= 0;
var y <= 0;
var z >= 0;

minimize Obj: x + 2*y + 3*z;

s.t. Con1:  x + y/z = 5;
s.t. Con2: -y + z/x = 7;
s.t. Con3:  z + x/y >=9;
