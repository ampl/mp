######################################
## dom_01.mod
## Test domain expressions
######################################

var x in -1..3 union 5..7;    ## AMPL computes upper bound 5
var y in 0..2 union interval [3,4];
var z in {0} union interval [5, 8.3];
var z1 in interval [2.2, 2.45] union interval [3.36, 6.8];
var z2 in union {i in 1..3} interval [i*3+0.1, i*3+2];

s.t. C1: x+y+z+z1+z2 == 11;

maximize O1: x;