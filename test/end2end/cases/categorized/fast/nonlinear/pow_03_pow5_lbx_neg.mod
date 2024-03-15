#########################################
## pow_03_pow5_lbx_neg.mod.
## Tests handling of pow(x) with lb(x)<0.
#########################################

var x >= -10, <= -2;
var y;

s.t. Pow01: y <= x^5;

s.t. Lin01: x+y <= -3;

maximize Obj01: 0.005*x+y;
