#########################################
## pow_01a_pow4_x_neg.mod.
## Tests handling of pow(x) with x<0.
#########################################

var x >= -10, <= 11;     ## Reduce range for COPT quadratics
var y;

s.t. Pow01: y >= x^4;

s.t. Lin01: x-y <= -3;

maximize Obj01: -y;
