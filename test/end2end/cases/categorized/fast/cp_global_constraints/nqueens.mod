
# -------------------------------------------------------------
# N-QUEENS (CSPLib prob054)
# nqueens.mod: place n queens on an n x n board so that none
# attacks another. Row[j] is the row of the queen in column j;
# three "alldiff" constraints: rows, and both diagonals.
# -------------------------------------------------------------

param n integer > 0 default 8;

set COLS := 1..n;

var Row {COLS} integer >= 1, <= n;

subj to RowAttacks:       alldiff {j in COLS} Row[j];
subj to UpDiagAttacks:    alldiff {j in COLS} (Row[j] + j);
subj to DownDiagAttacks:  alldiff {j in COLS} (Row[j] - j);

# Objective: a unique optimum that depends on the whole model
# (the weights make the board and its mirror image different)
minimize WeightedRows: sum {j in COLS} j * Row[j];
