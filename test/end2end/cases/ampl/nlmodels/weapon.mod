# weapon.mod	OLR2-MN-100-147
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Weapon assignment problem

# Ref.: D. M. Himmelblau, Applied Nonlinear Programming,
# McGraw-Hill Book Company, New York, 1972.  Problem 23.

# Number of (integer) variables:  100
# Number of constraints: 112 (147 after presolve adds more bounds)
# Objective nonseparable
# Objective nonconvex
# Linear constraints


param N integer > 0;
param M integer > 0;

set I := 1 .. N;
set J := 1 .. M;
set K;

param b {K} >= 0;
param c {I} <= 0;
param u {J} >= 0;
param a {I,J} >= 0, <= 1;

var x {I,J} integer >= 0;

minimize f:
   sum {j in J} u[j] * (prod {i in I} a[i,j]^x[i,j] - 1);

s.t. cc {j in K}:
    sum {i in I} x[i,j] >= b[j];

s.t. bb {i in I}:
    -sum {j in J} x[i,j] >= c[i];

data;
param N := 5;
param M := 20;
set K := 1 6 10 14 15 16 20;
param a  (tr):
           1        2        3       4       5     :=
      1    1       .84      .96      1       .92
      2   .95      .83      .95      1       .94
      3    1       .85      .96      1       .92
      4    1       .84      .96      1       .95
      5    1       .85      .96      1       .95
      6   .85      .81      .90      1       .98
      7   .90      .81      .92      1       .98
      8   .85      .82      .91      1        1
      9   .80      .80      .92      1        1
     10    1       .86      .95     .96      .90
     11    1        1       .99     .91      .95
     12    1       .98      .98     .92      .96
     13    1        1       .99     .91      .91
     14    1       .88      .98     .92      .98
     15    1       .87      .97     .98      .99
     16    1       .88      .98     .93      .99
     17    1       .85      .95      1        1
     18   .95      .84      .92      1        1
     19    1       .85      .93      1        1
     20    1       .85      .92      1        1    ;

param c :=
      1  -200
      2  -100
      3  -300
      4  -150
      5  -250 ;

param b :=
      1  30
      6  100
     10  40
     14  50
     15  70
     16  35
     20  10  ;

param u :=
      1  60     5  40     9  25    13  125    17  100
      2  50     6  60    10  150   14  200    18  100
      3  50     7  35    11  30    15  200    19  100
      4  75     8  30    12  45    16  130    20  150  ;
