###################################################
## Instance (11) from The3 SCIP Optimization Suite 8.0
###################################################

var x;
var y;

minimize Z:  exp( log(1000) + 1 + x*y );

s.t. CQ1: x^2 + y^2 <= 2;
