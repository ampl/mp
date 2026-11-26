# mip:basis after a MIP stopped at a limit (e.g., lim:mipsol=1 for
# Xpress): the fixed LP must still run and give a basis.  Xpress leaves
# the stopped MIP presolved, and its copy for the fixed LP refused the
# incumbent (?902 XPRSloadmipsol), so the report had no basis or duals.

param n := 60; param m := 10;
param w {i in 1..m, j in 1..n} := 1 + floor(Uniform(0, 100));
param c {i in 1..m} := 0.5 * sum {j in 1..n} w[i,j];
param p {j in 1..n} := sum {i in 1..m} w[i,j] / m + floor(Uniform(0, 20));

var x {1..n} binary;
var s >= 0, <= 1;

maximize v: sum {j in 1..n} p[j]*x[j] + s;

subject to cap {i in 1..m}: sum {j in 1..n} w[i,j]*x[j] + s <= c[i];
