# ljcluster.mod       OUR2-AY-V20-0
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Lennard-Jones Clusters in d dimensions

# Ref.: The MINPACK-2 Test Problem Collection, by B.M. Averick,
# R.G. Carter, J.J. More' and G.-L. Xue. Argonne National Laboratory,
# Mathematics and Computer Science Division, Preprint MCS-P153-0692,
# June 1992, pp. 41-45.

# Number of variables: 20;
# Number of constraints: 0;
# Objective partially separable
# Objective nonconvex

param d > 0 integer default 2;	# space dimension
set D := 1..d;
param N > 0 integer default 10;	# number of atoms
set I := {1..N};
set P := {i in I, j in 1..i-1};	# pairs of atoms

var x {i in I, D} default i;
var r {(i,j) in P} =		# distance separating atoms i and j
      sqrt( sum{k in D} (x[i,k] - x[j,k])^2 );

minimize energy:
      sum {(i,j) in P} (r[i,j]^-12 - 2*r[i,j]^-6);
