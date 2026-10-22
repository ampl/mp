
# -------------------------------------------------------------
# multi_sol_03.mod: multi_sol_02 with continuous overtime,
# so that the fixed LP of mip:basis has nonzero duals.
# With Xpress sol:poolmode=2 and mip:basis=1, the integer
# variables used to be fixed to the enumerator's last solution,
# giving duals inconsistent with the reported solution.
# -------------------------------------------------------------

set FAC := 1..4;
set CUST := 1..14;
param c {FAC, CUST};
param a {FAC, CUST};
param f {FAC};
param cap {FAC};

var x {FAC, CUST} binary;
var y {FAC} binary;
var over {FAC} >= 0, <= 6;

minimize Cost:
  sum {i in FAC, j in CUST} c[i,j]*x[i,j] + sum {i in FAC} f[i]*y[i]
  + sum {i in FAC} 2*over[i];

s.t. Assign {j in CUST}: sum {i in FAC} x[i,j] = 1;
s.t. Cap {i in FAC}: sum {j in CUST} a[i,j]*x[i,j] <= cap[i]*y[i] + over[i];

data;

param c (tr):
     1    2    3    4  :=
1    23   19   26   38
2    25   40   14   30
3    11   22   14   14
4    47   39   46   39
5    26   40   49   19
6    40   13   11   48
7    19   12   46   23
8    40   45   20   37
9    28   49   35   28
10   28   25   43   43
11   27   46   16   22
12   38   12   11   49
13   42   22   22   27
14   40   10   37   29
;

param a (tr):
     1    2    3    4  :=
1    15    9   13    5
2    15   11    6   16
3     6    6    5   10
4    17   17    8    9
5    17   10    5   18
6    11   11   13   14
7    13   10    8    6
8    19   14    5   15
9    17   12   12   15
10    6   13    7    5
11    6   16   13   19
12   16   19   14   14
13    7   10   19   14
14    9   12   16    7
;

param f := 1 22  2 55  3 40  4 31;
param cap := 1 52.2  2 51  3 43.2  4 50.1;

model;
