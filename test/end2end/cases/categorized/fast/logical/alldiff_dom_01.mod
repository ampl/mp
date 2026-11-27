# -------------------------------------------------------------
# alldiff_dom_01.mod
# alldiff over integer variables with 2100 values each: within
# cvt:uenc:maxdom (10000), the unary encoding is used and the
# model solves quickly.
# Optimum sum i*x[i] with x a permutation of 1..8 in decreasing
# order: 120.
# -------------------------------------------------------------
param n := 8;

var x {1..n} integer >= 1 <= 2100;

minimize o: sum {i in 1..n} i*x[i];

s.t. c: alldiff {i in 1..n} x[i];
