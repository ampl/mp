# levy3.mod	OBR2-AN-2-4
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: C. Jansson and O. Knueppel, "A Global Minimization Method:
# the Multi-Dimensional Case", Technische Informatik III,
# TU Hamburg-Hamburg, Jan. 1992, pp. 71-72 (problem "whs4").

# Levy function No. 3

# Number of variables:  2
# Number pf constraints:  4
# Objective separable
# Objective nonconvex
# Simple bound constraints

# Note: there is a mistake in the values of global minima in the Ref.

# There are lots of local minima in the domain.

set I := 1 .. 2;
var x{I} <= 10, >= -10;

minimize Flevy:
  - prod{i in I} sum {j in 1..5} (j*cos((j + 1)*x[i] + j));

data;

var x  :=
    1  4.976
    2  4.976;
