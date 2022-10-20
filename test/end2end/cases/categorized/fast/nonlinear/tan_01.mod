#########################################
## cos_01.mod
## Check an intermediate point of tangens
#########################################

param pi := 4*atan(1);  ## 3.14159265358979;

var x >= -4*pi, <= 1*pi;
var y;

s.t. Tan01: y >= tan(x);

s.t. Lin01: y <= -8.235*(x+0.9333*pi);

s.t. Lin02: y >= -0.5235*(x+2.16*pi);

maximize Obj01: 0.5*x-y;
