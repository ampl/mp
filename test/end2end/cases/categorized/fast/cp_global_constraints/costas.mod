
# -------------------------------------------------------------
# COSTAS ARRAY
# costas.mod: find a permutation P of 1..n such that, for every
# distance d, the differences P[i+d] - P[i] are all different.
# (Equivalently, no two vectors between pairs of dots coincide.)
# Uses a family of "alldiff" constraints over expressions.
# -------------------------------------------------------------

param n integer > 1 default 8;

set POS := 1..n;

var P {POS} integer >= 1, <= n;

subj to Perm:  alldiff {i in POS} P[i];

subj to DistinctDiffs {d in 1..n-2}:
   alldiff {i in 1..n-d} (P[i+d] - P[i]);

# Symmetry breaking (mirror image)
subj to Mirror:  P[1] < P[n];

# Objective: the largest possible P[1]
maximize First: P[1];
