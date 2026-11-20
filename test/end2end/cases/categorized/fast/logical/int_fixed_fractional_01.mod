# -------------------------------------------------------------
# int_fixed_fractional_01.mod
# An integer variable fixed at a fractional value is infeasible.
# Wrong before (without the AMPL presolve): the fixed variable was
# made continuous and the model "solved" with x = 0.5.
# -------------------------------------------------------------

var x integer >= -1 <= 2;

minimize o: x;

s.t. c1: 2*x == 1;
