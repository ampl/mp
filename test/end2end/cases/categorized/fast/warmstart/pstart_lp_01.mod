# Primal start values without dual values (no constraint has a value).
# Optimum 17: x[2] = 5, x[4] = 7.
# The suffix is declared here, so that it reads 0 when not reported.
suffix test_ini_pri OUT;

set J := 1..4;
param c {j in J} := if j = 1 then 3 else if j = 2 then 2
                      else if j = 3 then 4 else 1;

var x {j in J} >= 0 <= 10 := if j = 2 then 4 else if j = 4 then 8 else 0;

minimize cost: sum {j in J} c[j] * x[j];

s.t. demand: sum {j in J} x[j] >= 12;
s.t. mix:    x[1] + x[2] >= 5;
s.t. cap:    x[3] + x[4] <= 9;
