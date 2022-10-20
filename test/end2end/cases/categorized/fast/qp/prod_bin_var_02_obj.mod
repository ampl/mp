#########################
# prod_bin_var_02_obj.mod
# Product with a binary in the objective
#########################

var b: binary;
var y >=-15, <=6;

minimize Prod01: b*(y+3)+4;

s.t. Sum01: b+y >= -5;
