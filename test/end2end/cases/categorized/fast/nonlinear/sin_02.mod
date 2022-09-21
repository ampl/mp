#########################################
## sin_02.mod
## Check the left end of a sinus curve
#########################################

param pi := 4*atan(1);  ## 3.14159265358979;

var x >= 1.5*pi, <= 15*pi;
var y >=-1, <=1;

s.t. Sin01: y <= sin(x);

maximize SumXY: -x+y;
