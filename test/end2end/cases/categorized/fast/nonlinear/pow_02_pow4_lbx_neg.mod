#########################################
## pow_02_pow4_lbx_neg.mod.
## Tests handling of pow(x) with lb(x)<0.
#########################################

var x >= -10, <= 20;
var y;

s.t. Pow01: y >= x^4;

s.t. Lin01: -x+y >= 3;

maximize Obj01: 0.005*x-y;
