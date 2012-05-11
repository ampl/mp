# schwefel.mod	OBR2-AN-5-10
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: C. Jansson and O. Knueppel, "A Global Minimization Method:
# the Multi-Dimensional Case", Technische Informatik III,
# TU Hamburg-Hamburg, Jan. 1992, p. 103 (problem "whs25").

# Schwefel function

# Number of variables:  5
# Number of constraints:  10
# Objective separable convex
# Simple bound constraints

# The global minimum is Fs = 0, x = (0,0,0,0,0).

set I := {1..5};

var x{I} <= 0.4, >= -0.5, := .4;

minimize Fs:
     sum {i in I} x[i]^10;
