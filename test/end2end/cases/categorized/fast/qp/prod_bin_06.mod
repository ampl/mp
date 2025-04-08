##########################################
## Test products of binary variables, in particular reformulations into logicals
## of 2-term products.
## prod_bin_06.mod
##########################################

var b {1..5} binary;

minimize Obj: b[1] - b[2] + b[3] - b[4] + b[5];

s.t. C1: b[1]*((b[2]-1)*(b[3]-1)-1)*b[4] == -1;

s.t. C2: b[2]*b[3] + (1-b[3])*b[5] + b[2]*(1-b[5])  ## To use QP2Pass #153
          + 6*b[3] + 7*b[2] >= 5;