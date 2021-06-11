# hs6.mod		QQR2-RN-2-1
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 29.

# Number of variables: 2
# Number of constraints: 1
# Objective quadratic
# Nonlinear constraint

var x{1..2};

minimize f:
         (1 - x[1])^2;

s.t. Constr:
         10 * (x[2] - x[1]^2) = 0;


data;
var x:=
    1 -1.2
    2  1 ;
