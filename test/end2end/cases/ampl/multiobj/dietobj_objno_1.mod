set MINREQ;   # nutrients with minimum requirements
set MAXREQ;   # nutrients with maximum requirements

set NUTR := MINREQ union MAXREQ;   # nutrients
set FOOD;                          # foods
set STORE;                         # stores

param cost {STORE,FOOD} > 0;
param f_min {FOOD} >= 0;
param f_max {j in FOOD} >= f_min[j];

param n_min {MINREQ} >= 0;
param n_max {MAXREQ} >= 0;

param amt {NUTR,FOOD} >= 0;

var Buy {j in FOOD} >= f_min[j], <= f_max[j];

minimize total_cost {s in STORE}:
   sum {j in FOOD} cost[s,j] * Buy[j];

minimize total_number:  sum {j in FOOD} Buy[j];

subject to diet_min {i in MINREQ}:
   sum {j in FOOD} amt[i,j] * Buy[j] >= n_min[i];

subject to diet_max {i in MAXREQ}:
   sum {j in FOOD} amt[i,j] * Buy[j] <= n_max[i];
