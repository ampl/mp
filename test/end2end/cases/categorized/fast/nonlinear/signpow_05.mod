## Test abs(C*x) * x

var x >= -5 <= 3;
var y;

param ObjTerm1Coef default 0.1;
param C default -5;
param ConXOffset default -2;
param ConXTermPow default 3;

minimize SPDiff: ObjTerm1Coef * abs(C*x) * x + exp(x) - y;

s.t. SP2: y <= abs(x + ConXOffset) * (x + ConXOffset)^ConXTermPow;
