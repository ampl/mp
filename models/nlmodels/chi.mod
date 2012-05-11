# chi.mod OBR2-AN-2-4
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: C. Jansson and O. Knueppel, "A Global Minimization Method:
# the Multi-Dimensional Case", Technische Informatik III,
# TU Hamburg-Hamburg, Jan. 1992, p. 67 (problem "Chi").

# Chichinadze function

# Number of variables: 2
# Number of constraints: 4
# Objective separable
# Objective nonconvex
# Simple bound constraints

# The global minimum is Fchi = -43.31586207 at x = (5.90133,0.5).

param pi := 4*atan(1);
var x{1..2};

minimize Fchi:
	  x[1]^2 - 12*x[1] + 11 + 10*cos(pi*x[1]/2)
	+ 8*sin(pi*5*x[1]) - exp(-(x[2] - .5)^2/2)/sqrt(5);

s.t. Box1:
     -30 <= x[1] <= 30;

s.t. Box2:
     -10 <= x[2] <= 10;
