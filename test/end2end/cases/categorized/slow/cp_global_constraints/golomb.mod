
# -------------------------------------------------------------
# GOLOMB RULER
# golomb.mod: place n marks on a ruler so that all distances
# between pairs of marks are different, minimizing the length
# of the ruler. Uses "alldiff" over the differences.
#
# Known optimal lengths (OEIS A003022):
#   n    : 4  5  6  7  8  9  10 11 12  13
#   len  : 6  11 17 25 34 44 55 72 85  106
# -------------------------------------------------------------

param n integer > 1;

set MARKS := 1..n;

# n^2 is a (loose) upper bound on the optimal length
var Mark {MARKS} integer >= 0, <= n^2;

minimize Length: Mark[n];

subj to FirstMarkAtZero: Mark[1] = 0;

subj to Increasing {i in 1..n-1}: Mark[i] + 1 <= Mark[i+1];

subj to DistinctDistances:
   alldiff {i in MARKS, j in MARKS: i < j} (Mark[j] - Mark[i]);

# Symmetry breaking: a ruler and its mirror image are the same
subj to Mirror: Mark[2] - Mark[1] < Mark[n] - Mark[n-1];
