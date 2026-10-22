# -------------------------------------------------------------
# IIS LogicalCon 08
# infeas_lcon_08.mod
# Every constraint uses 'if ... then', whose else-branch is
# the constant 0. The infeasibility is in period 3 only.
# Constraints of other periods must not be in the IIS:
# MP's shared fixed variable for the constant used to pass
# its IIS status to every constraint using it.
# -------------------------------------------------------------

set U := 1..2;
set T := 1..4;
param demand {t in T} := if t == 3 then 18 else 8;

var P {U, T} >= 0, <= 10;
var F {U, T} >= 0;

minimize Cost: sum {u in U, t in T} F[u,t];

s.t. Fuel {u in U, t in T}: F[u,t] >= 2*P[u,t] + (if P[u,t] >= 1 then 5);
s.t. Reserve {t in T}: sum {u in U} (if P[u,t] >= 1 then 10 - P[u,t]) >= 3;
s.t. Balance {t in T}: sum {u in U} P[u,t] == demand[t];
