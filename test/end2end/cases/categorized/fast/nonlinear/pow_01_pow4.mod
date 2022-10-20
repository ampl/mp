#########################################
## pow_01_pow4.mod
#########################################

var x >= -10, <= 11;     ## Reduce range for COPT quadratics
var y;

s.t. Pow01: y >= x^4;

s.t. Lin01: x+y >= 3;

maximize Obj01: 0.005*x-y;
