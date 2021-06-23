# gold.mod	OBR2-AN-2-4
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: C. Jansson and O. Knueppel, "A Global Minimization Method:
# the Multi-Dimensional Case", Technische Informatik III,
# TU Hamburg-Hamburg, Jan. 1992, p. 37 (problem "GP").

# Goldstein-Price function

# Number of variables:  2
# Number of constraints:  4
# Objective nonseparable
# Objective nonconvex
# Simple bound constraints


# The global minimum is Fgold = 3, x = (0,-1).

var x{1..2} <= 2, >= -2;

minimize Fgold:
   (1 + (x[1] + x[2] + 1)^2
	* (19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2))
 * (30 + (2*x[1] - 3*x[2])^2
	* (18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2));

data;

var x  :=
    1   0.5
    2  -0.5;
