# steenbre.mod	ONR2-MY-540-666
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# A nonconvex multi-commodity network problem

# Ref.: P. A. Steenbrink,Optimization of Transport Network,
# Wiley, 1974, p.120.

# The CUTE description of the same problem can be found in STEENBRE.SIF .
# The problem has been slightly perturbed by making cmin = 0.01 instead of 0,
# and by the addition of sqrt_offset to the sqrt terms.

# Number of variables: 540
# Number of constraints: 126
# Objective partially separable nonconvex
# Linear constraints

set CITIES := {1..9};
set DEST := {1..6};
set TRIPS within (DEST cross DEST);
set ROADS within (CITIES cross CITIES);
param cost {(i,j) in ROADS} >= 0;
check {(i,j) in ROADS: i < j}: cost[i,j] = cost[j,i];
param cmin {ROADS} > 0;
param TZERO >= 0;
param CCR >= 0;
param alpha {ROADS} >= 0;
param tr_matr {TRIPS};
param sqrt_offset default .01;	# from CUTE's STEENBRE.SIF

node Balance {(k,l) in TRIPS, j in CITIES}:
     net_in =
              if l = j  then tr_matr[k,l]
              else if k = j  then -tr_matr[k,l]
              else 0;

var capacity {(i,j) in ROADS} >= cmin[i,j];
arc flow {(k,l) in TRIPS, (i,j) in ROADS} >= 0,
    from Balance[k,l,i], to Balance[k,l,j];
var total_flow {(i,j) in ROADS} = sum {(k,l) in TRIPS} flow[k,l,i,j];

minimize Total_Cost:
         sum {(i,j) in ROADS} (
		  cost[i,j] * alpha[i,j] *
			sqrt(capacity[i,j] - cmin[i,j] + sqrt_offset)
		+ total_flow[i,j] * cost[i,j] *
			(TZERO + CCR * (total_flow[i,j] / capacity[i,j])^2)
		);



data;

param alpha	default .01;
var capacity	default .1 :=	[2,4] 2000.1	[4,2] 2000.1;
var flow	default .1;
set ROADS:   1   2   3   4   5   6   7   8   9   :=
      1      -   +   +   -   -   -   +   -   -
      2      +   -   -   +   -   -   +   +   -
      3      +   -   -   -   +   -   +   +   -
      4      -   +   -   -   -   +   -   +   +
      5      -   -   +   -   -   +   -   +   +
      6      -   -   -   +   +   -   -   -   +
      7      +   +   +   -   -   -   -   +   -
      8      -   +   +   +   +   -   +   -   +
      9      -   -   -   +   +   +   -   +   -    ;

param:  TRIPS:          tr_matr      :=
        1  6            10000
        2  3            2000
        2  4            2000
        2  5            1000
        3  2            200
        3  4            1000
        3  5            2000
        4  2            200
        4  3            100
        4  5            1000
        5  2            100
        5  3            200
        5  4            100
        6  1            1000   ;



param   cost:
         1     2    3    4    5    6    7    8    9           :=
   1     .    35   40    .    .    .   30    .    .
   2    35     .    .   100   .    .   15    55   .
   3    40     .    .    .   100   .   25    60   .
   4     .    100   .    .    .   35    .    55   15
   5     .     .   100   .    .   40    .    60   25
   6     .     .    .   35   40    .    .     .   30
   7    30    15   25    .    .    .    .    50   .
   8     .    55   60   55   60    .   50     .   50
   9     .     .    .   15   25   30    .    50   .   ;

param TZERO := 0.01;
param CCR := 0.01;

param cmin default 0.01 :=
	[2,4]  2000
	[4,2]  2000;
