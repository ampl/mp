
# -------------------------------------------------------------
# LANGFORD'S NUMBER PROBLEM (CSPLib prob024), L(2,n)
# langford.mod: arrange two copies of each of the numbers 1..n in
# a sequence of length 2n such that the two copies of k have
# exactly k numbers between them. A solution exists iff n mod 4
# is 0 or 3.
# Uses "alldiff" over all 2n positions.
# -------------------------------------------------------------

param n integer > 0 default 7;

set NUMS := 1..n;

# Pos[k,c] is the position of copy c of the number k
var Pos {k in NUMS, c in 1..2} integer >= 1, <= 2*n;

subj to Spacing {k in NUMS}:  Pos[k,2] = Pos[k,1] + k + 1;

subj to AllPositionsDiffer:   alldiff {k in NUMS, c in 1..2} Pos[k,c];

# Symmetry breaking: a sequence and its reversal are both solutions
subj to Reversal:  Pos[1,1] <= n;

# Objective: the earliest possible position of the first "1"
minimize FirstOne: Pos[1,1];
