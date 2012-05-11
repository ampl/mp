# AMPL model for an MPS file with 'MARKER' lines indicating integer
# variables: like mps.mod, this is short, but it may
# change the row order (which could affect pivot choice
# by simplex-based solvers).

# Use the awk script "m2ai" to turn an MPS file into suitable data.

set Aij dimen 2;		#constraint matrix indices
set I1;				# to allow empty rows
set J := setof{(i,j) in Aij} j;	#columns
set K within J;			#integer variables
param A{Aij};			#constraint matrix nonzeros

param b{I1} default 0;		#right-hand side
param db{I1};			#for ranges

set ctypes := {'N', 'L', 'E', 'G', 'LR', 'GR'};

param ctype{I1} symbolic within ctypes;

param lb{J} default 0;
param ub{J} default Infinity;

var x{j in J diff K} >= if lb[j] <= -1.7e38 then -Infinity else lb[j]  <= ub[j];

var ix{j in K} integer >= if lb[j] <= -1.7e38 then -Infinity else lb[j]
		  <= ub[j];

set zork := setof{i in I1} (i,ctype[i]);

var xx{j in J} = if j in K then ix[j] else x[j];

minimize Obj{(i,'N') in zork}:  sum{(i,j) in Aij} A[i,j]*xx[j];

E{(i,'E') in zork}:  sum{(i,j) in Aij} A[i,j]*xx[j] == b[i];
L{(i,'L') in zork}:  sum{(i,j) in Aij} A[i,j]*xx[j] <= b[i];
G{(i,'G') in zork}:  sum{(i,j) in Aij} A[i,j]*xx[j] >= b[i];
LR{(i,'LR') in zork}:  b[i] - db[i] <= sum{(i,j) in Aij} A[i,j]*xx[j] <= b[i];
GR{(i,'GR') in zork}:  b[i] + db[i] >= sum{(i,j) in Aij} A[i,j]*xx[j] >= b[i];

option linelim 1;
