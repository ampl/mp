#########################################
## tanh_01.mod
#########################################

param pi := 4*atan(1);  ## 3.14159265358979;

var x >= -5, <= 3;      ## HiGHS cannot handle ample bounds robustly
var y;

s.t. Sin01: y == tanh(x);

s.t. Lin01: -x-y <= 1;

maximize SumKXY: -13*x-y;
