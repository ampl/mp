###############################################
## Test repeated redefinition of a functional constraint
###############################################

var x binary;
var y>=-15 <=18;

maximize Obj: x+y;

s.t. C1: x - y < 19;

s.t. C2: exp(x*y) - x <= 8;