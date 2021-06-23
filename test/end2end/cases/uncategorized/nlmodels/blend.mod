# blend.mod	LOR2-MN-24-38
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Blending multicomponent mixtures

# Ref.: D. M. Himmelblau, Applied Nonlinear Programming,
# McGraw-Hill Book Company, New York, 1972.  Problem 20.

# Number of variables:  24
# Number of constraints:  38
# Objective nonseparable
# Objective nonconvex
# Nonlinear constraints


param ncomp > 0  integer;
param neq > 0 integer;
set I := {1..ncomp};
set J := {1..neq};
param a {i in I} > 0;
param b {i in I} > 0;
param c {i in J} > 0;
param d {i in J} > 0;
param f1 > 0;
param f2 > 0;
param f3 > 0;
param f4 > 0;
param f := f1*f2*f3/f4;


var x {i in I} >= 0;

minimize cost:
         sum {i in I} a[i] * x[i];

s.t. eq {i in J} :
     x[i+neq] / (b[i+neq] * sum {j in J} x[j+neq] / b[j+neq])
     ==
     c[i] * x[i] / (40 * b[i] * sum {j in J} x[j] / b[j]);

s.t. simplex:
     sum {i in I} x[i] == 1;

s.t. nl2:
     sum {i in J} x[i] / d[i] + f * sum {i in J} x[i+neq] / b[i+neq] == 1.671;

## Himmelblau (p. 420) gives two other sets of nonlinear inequality
## constraints, here omitted, which are apparently stated incorrectly:
## they imply x[i] == 0 for i in {1,2,3,7,8,9,13,14,15,19,20,21}.

data;

param ncomp := 24;
param neq := 12;
param :          a        b      :=
        1      0.0693   44.094
        2      0.0577   58.12
        3      0.05     58.12
        4      0.20    137.4
        5      0.26    120.9
        6      0.55    170.9
        7      0.06     62.501
        8      0.10     84.94
        9      0.12    133.425
       10      0.18     82.507
       11      0.10     46.07
       12      0.09     60.097
       13      0.0693   44.094
       14      0.0577   58.12
       15      0.05     58.12
       16      0.20    137.4
       17      0.26    120.9
       18      0.55    170.9
       19      0.06     62.501
       20      0.10     84.94
       21      0.12    133.425
       22      0.18     82.507
       23      0.10     46.07
       24      0.09     60.097 ;                                     ;

param:            c       d      :=
        1     123.7    31.224
        2      31.7    36.12
        3      45.7    34.784
        4      14.7    92.7
        5      84.7    82.7
        6      27.7    91.6
        7      49.7    58.708
        8       7.1    82.7
        9       2.1    80.8
       10      17.7    64.517
       11       0.85   49.4
       12       0.64   49.1  ;


param f1 := 0.7302;
param f2 := 530;
param f3 := 14.7;
param f4 := 40;

var:       x        :=
     1    0.04
     2    0.04
     3    0.04
     4    0.04
     5    0.04
     6    0.04
     7    0.04
     8    0.04
     9    0.04
    10    0.04
    11    0.04
    12    0.04
    13    0.04
    14    0.04
    15    0.04
    16    0.04
    17    0.04
    18    0.04
    19    0.04
    20    0.04
    21    0.04
    22    0.04
    23    0.04
    24    0.04;
