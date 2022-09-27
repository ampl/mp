#########################################
## pow_05_pow_minus0.35.mod
#########################################

var x >= 0.4, <= 1e2;   ## HiGHS does not like large bound
var y;

s.t. Pow01: y >= x^-.35;

s.t. Lin01: x+y >= 3;

minimize Obj01: x+2*y;
