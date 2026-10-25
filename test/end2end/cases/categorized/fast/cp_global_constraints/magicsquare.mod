
# -------------------------------------------------------------
# MAGIC SQUARE
# magicsquare.mod: fill an n x n square with the numbers 1..n^2,
# each used once, so that all rows, columns and both diagonals
# sum to n*(n^2+1)/2.
# Combines "alldiff" with linear constraints.
# -------------------------------------------------------------

param n integer > 2 default 4;

set R := 1..n;

param magic := n * (n^2 + 1) / 2;

var X {R, R} integer >= 1, <= n^2;

subj to AllDiff:  alldiff {i in R, j in R} X[i,j];

subj to Rows {i in R}:  sum {j in R} X[i,j] = magic;
subj to Cols {j in R}:  sum {i in R} X[i,j] = magic;
subj to Diag1:          sum {i in R} X[i,i] = magic;
subj to Diag2:          sum {i in R} X[i,n+1-i] = magic;

# Symmetry breaking
subj to Sym1:  X[1,1] < X[1,n];
subj to Sym2:  X[1,n] < X[n,1];
subj to Sym3:  X[1,1] < X[n,n];

# Objective: largest possible value in the top-left corner (the
# smallest of the four corners, by the symmetry breaking above)
maximize Corner: X[1,1];
