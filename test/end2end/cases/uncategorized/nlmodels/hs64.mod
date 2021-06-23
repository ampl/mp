# hs64.mod	OOR2-MN-3-4
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 86.

# Number of variables: 3
# Number of constraints: 4
# Objective separable
# Objective nonconvex
# Nonlinear constraints

var x{1..3} >= 1E-5, := 1;

minimize Obj:
         5 * x[1] + 50000 / x[1] + 20 * x[2] + 72000 / x[2] + 10 * x[3] +
         144000 / x[3];

s.t. Constr:
     1 >= 4/x[1] + 32/x[2] + 120/x[3];
