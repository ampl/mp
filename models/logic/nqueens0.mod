
param n;

set ROWS := {1..n};
set COLUMNS := {1..n};

var X{ROWS, COLUMNS} binary; 
   # X[i,j] is one if there is a queen at (i,j); else zero

#maximize max_queens: sum {i in ROWS, j in COLUMNS} X[i,j];

subject to column_attacks {j in COLUMNS}:
	sum {i in ROWS} X[i,j] = 1;

subject to row_attacks {i in ROWS}:
	sum {j in COLUMNS} X[i,j] = 1;

subject to diagonal1_attacks {k in 2..2*n}:
	sum {i in ROWS, j in COLUMNS: i+j=k} X[i,j] <= 1;

subject to diagonal2_attacks {k in -(n-1)..(n-1)}:
	sum {i in ROWS, j in COLUMNS: i-j=k} X[i,j] <= 1;
