#########################################
## cos_01.mod
## Check an intermediate point of cosine
#########################################

param pi := 4*atan(1);  ## 3.14159265358979;

var x >= 40*pi, <= 81*pi;
var y >=-1, <=1;

s.t. Cos01: y == cos(x);

s.t. Lin01: y <= 5*(x-61.333*pi);

maximize Obj01: 3*y-x;
