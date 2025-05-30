param n_pieces default 4;
param n_towers default 3;
set PIECES := 1..n_pieces;
set TOWERS := 1..n_towers;
set POSITIONS := 1..n_pieces;
set TURNS := 1..(2^n_towers + 3); # just in case
param weight {PIECES} default 1;

var x{PIECES, POSITIONS, TOWERS, TURNS} binary;

init{i in PIECES}: x[i,i,1,1] = 1;

one_piece_one_tower {i in PIECES, t in TURNS}:
	sum{j in POSITIONS, k in TOWERS} x[i,j,k,t] = 1;
stg_above2
		{i in PIECES, j in POSITIONS, k in TOWERS, t in TURNS: t < 2^n_towers}:
	sum{ii in PIECES, jj in POSITIONS: i != ii and jj > j} x[ii,jj,k,t] >= 1 ==> x[i,j,k,t] = x[i,j,k,t+1];
final{i in PIECES}: x[i,1,n_towers,2^n_towers + 3] = 1;
