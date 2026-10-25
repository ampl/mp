
# -------------------------------------------------------------
# ALL-INTERVAL SERIES (CSPLib prob007)
# allinterval.mod: find a permutation x of 0..n-1 such that the
# absolute differences between consecutive elements are also all
# different (so they are a permutation of 1..n-1).
# Uses "alldiff" over expressions containing abs().
# -------------------------------------------------------------

param n integer > 1 default 9;

set POS := 1..n;

var X {POS} integer >= 0, <= n-1;

subj to Perm:  alldiff {i in POS} X[i];

subj to Intervals:
   alldiff {i in 1..n-1} abs (X[i+1] - X[i]);

# Objective: a unique optimum that depends on the whole model
maximize WeightedPositions: sum {i in POS} i * X[i];
