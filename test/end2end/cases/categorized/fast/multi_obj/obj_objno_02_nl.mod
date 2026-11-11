#######################################
## obj_objno_02_nl.mod
## objno selects the 2nd of two nonlinear objectives
## with different minimizers and constant terms
#######################################

var x >= 2.5 <= 3.5 := 2.6;
var y >= 0.5 <= 1.5 := 0.6;

minimize f1: exp(x-2.7) - x + (y-1.3)^4 + 7;   ## optimum 5.3 at (2.7, 1.3)
minimize f2: exp(x-3) - x + (y-1)^4 + 3;       ## optimum 1 at (3, 1)

s.t. C: x + y <= 5;
