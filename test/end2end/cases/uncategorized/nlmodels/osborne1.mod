# osborne1.mod OBR2-MN-5-10
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Exponential data fitting

# Original reference:
# M. R. Osborne, "Some Aspects of Nonlinear Least Squares Calculations".
# In Numerical Methods for Nonlinear Optimization, edited by F. A.
# Lootsma, Academic Press, New York, 1972.

# This also appears as problem "Exponential Data Fitting 1" in
# "The MINPACK-2 Test Problem Collection", by B.M. Averick,
# R.G. Carter, J.J. More' and G.-L. Xue; Argonne National Laboratory,
# Mathematics and Computer Science Division, Preprint MCS-P153-0692,
# June 1992, p. 26-27.
# The bounds are from this latter reference.

# Number of variables: 5
# Number of constraints: 10
# Objective nonseparable
# Objective nonconvex
# Simple bound constraints

# Global minimum ssq = 5.464894697e-05 at
# x = (0.37541, 1.93585, -1.46469, 0.0128675, 0.0221227).

param N > 0 integer;
param M > 0 integer;
set I := 1 .. N;
set J := 1 .. M;
param y {j in J};
param t {j in J} := 10 * (j - 1);

var x {i in I} >= -10, <= 10;

minimize ssq:
   sum {j in J}
   (y[j] - (x[1] + x[2]*exp(-t[j]*x[4]) + x[3]*exp(-t[j]*x[5])))^2;


data;
param N := 5;
param M := 33;
param y :=
    1 .844    7 .881   13 .685   19 .538   25 .457   31 .414
    2 .908    8 .85    14 .658   20 .522   26 .448   32 .411
    3 .932    9 .818   15 .628   21 .506   27 .438   33 .406
    4 .936   10 .784   16 .603   22 .49    28 .431
    5 .925   11 .751   17 .58    23 .478   29 .424
    6 .908   12 .718   18 .558   24 .467   30 .42
    ;

var x :=	# initial guess
	1   .5
	2  1.5
	3 -1
	4   .01
	5   .02
	;
