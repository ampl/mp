
# -------------------------------------------------------------
# MAGIC SEQUENCE (CSPLib prob019)
# magicseq.mod: find a sequence x[0..n-1] such that x[i] is the
# number of occurrences of the value i in the sequence.
# For n >= 7 the solution is unique: x[0] = n-4, x[1] = 2,
# x[2] = 1, x[n-4] = 1, all other elements are 0.
# Uses the "numberof" operator with a constant, and "count".
# -------------------------------------------------------------

param n integer > 3 default 10;

set IDX := 0..n-1;

var X {IDX} integer >= 0, <= n-1;

subj to Magic {i in IDX}:
   X[i] = numberof i in ({j in IDX} X[j]);

# Redundant constraints: sum of x[i] is n
subj to SumIsN:  sum {i in IDX} X[i] = n;

# Same counting, with "count": the number of j's with X[j] = 0
# equals X[0].
subj to ZerosCount:  X[0] = count {j in IDX} (X[j] = 0);

maximize FirstElement: X[0];
