# branin.mod OBR2-AN-2-4
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: C. Jansson and O. Knueppel, "A Global Minimization Method:
# the Multi-Dimensional Case", Technische Informatik III,
# TU Hamburg-Hamburg, Jan. 1992, p. 35 (problem "BR").

# Number of variables: 2
# Number of constraints:  4
# Objective nonseparable, nonconvex
# Simple bound constraints


# There are 3 global minima with Fbr = 0.397887: (-pi,12.275),
# (pi,2.275), (3*pi,2.475).

param pi := 4*atan(1);
var x{1..2};

minimize Fbr:
	  (x[2] - 5.1*x[1]^2/(4*pi*pi) + 5*x[1]/pi - 6)^2
	+ 10*(1 - 1/(8*pi))*cos(x[1]) + 10;

s.t. Box1:
     -5 <= x[1] <= 10;

s.t. Box2:
     0 <= x[2] <= 15;
