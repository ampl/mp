# kowalik.mod	OBR2-AN-4-8
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Kowalik & Osborne problem, originally stated in

# J. R. Kowalik and M. R. Osborne, Methods for Unconstrained Minimization,
# American Elsevier, New York, 1968.

# but with modified b[i] values as specified in

# C. Jansson and O. Knueppel, "A Global Minimization Method: the Multi-
# Dimensional Case", Technische Informatik III, TU Hamburg-Hamburg, Jan. 1992.

# where the problem is called just "Kowalik".  The original formulation
# explicitly gives the b[i] values to three significant figures.


# Number of variables:  4
# Number of constraints:  8
# Objective nonconvex, nonseparable
# Simple bound constraints

# The global minimum is Fko = 3.07486*10e-4,
# x = (0.192833, 0.190836, 0.123117, 0.135766).

set I := {1..11};
param a{I};
param c{I};
param b{i in I} := 1/c[i];

var x{1..4} <= 0.42, >= 0, := .42;

minimize Fko:
         sum {i in I} (a[i] - x[1] * (b[i]^2 + b[i]*x[2]) /
                                     (b[i]^2 + b[i]*x[3] + x[4]))^2;

data;

param:      a        c  :=
      1   0.1957     0.25
      2   0.1947     0.5
      3   0.1735     1
      4   0.16       2
      5   0.0844     4
      6   0.0627     6
      7   0.0456     8
      8   0.0342    10
      9   0.0323    12
     10   0.0235    14
     11   0.0246    16  ;
