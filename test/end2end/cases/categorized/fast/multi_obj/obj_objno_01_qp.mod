#######################################
## obj_objno_01_qp.mod
## objno selects the 2nd of two quadratic objectives
## with different minimizers and constant terms
#######################################

var x >= 0 <= 4 := 1;
var y >= 0 <= 4 := 1;

minimize f1: (x-1)^2 + (y-2)^2 + 7;    ## optimum 7 at (1, 2)
minimize f2: (x-3)^2 + y^2 + 1;        ## optimum 1 at (3, 0)

s.t. C: x + y <= 5;
