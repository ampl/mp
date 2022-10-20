#########################
# prod_0_obj.mod
# Product in the objective
#########################

var x >=-3, <=7;
var y >=-1, <=6;

minimize Prod01: (x-8)*(y+3)+4;

s.t. Sum01: (x-8) - 2*(y+3) >= -10;

s.t. Sum02: 2*(x-8) - (y+3) >= -13;
