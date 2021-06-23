# s383.mod	OLR2-AN-14-29
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: K. Schittkowski, More Test Examples for Nonlinear Programming Codes.
# Lecture Notes in Economics and Mathematical Systems, v. 282,
# Springer-Verlag, New York, 1987, p. 202.

# Number of variables:  14
# Number of constraints:  29
# Objective separable convex
# Linear constraints

set I := 1 .. 14;

param a{I};
param c{I};

var x{I} >= 0 := 0.01;

minimize Obj:
         sum {i in I} a[i] / x[i];

s.t. G1:
     sum {i in I} c[i] * x[i] == 1;
s.t. B1 {i in {1..5}}:
     x[i] <= 0.04;
s.t. B2 {i in {6..14}}:
     x[i] <= 0.03;

data;

param a :=
      1    12842.275       8    1267
      2      634.25        9     760.05
      3      634.25       10     633.25
      4      634.125      11    1266.25
      5     1268          12     632.875
      6      633.875      13     394.46
      7      633.75       14     940.838   ;

param c :=
      1     5.47934        8     3.90896
      2      .83234        9     2.74284
      3      .94749       10     2.60541
      4     1.11082       11     5.96184
      5     2.64824       12     3.29522
      6     1.55868       13     1.83517
      7     1.73215       14     2.81372   ;
