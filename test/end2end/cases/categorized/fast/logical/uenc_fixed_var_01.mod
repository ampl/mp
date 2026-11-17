# -------------------------------------------------------------
# uenc_fixed_var_01.mod
# Seven guests, three tables, at most three guests at a table,
# every table needs two acquainted pairs (count >= 2); the
# acquaintance graph is a tree, so this is infeasible. Guest 1
# is fixed at table 1 by a separate constraint after the
# comparisons tables[1] == 2 and == 3 were created: their flags
# were left free and counted. Wrong before: "optimal".
# Reduced from Hakan Kjellerstrand's wedding_optimal_chart.mod
# (MIT license, see the hakank folder of the uncategorized cases).
# -------------------------------------------------------------

param n; param a; param b; param m; param C{1..m, 1..m};
var tables{1..m} >= 1 <= n integer;
minimize o: 0;
s.t. c2a{i in 1..n}: count{j in 1..m, k in 1..m: j < k} (C[j,k] > 0 and tables[j] = i and tables[k] = i) >= b;
s.t. c2b{i in 1..n}: sum{j in 1..m} (if tables[j] = i then 1) <= a;
s.t. c3: tables[1] = 1;
data;
param m := 7; param n := 3; param a := 3; param b := 2;
param C: 1 2 3 4 5 6 7 :=
1 0 0 0 0 1 1 1
2 0 0 0 0 0 0 1
3 0 0 0 0 1 0 0
4 0 0 0 0 0 1 0
5 1 0 1 0 0 0 0
6 1 0 0 1 0 0 0
7 1 1 0 0 0 0 0;
