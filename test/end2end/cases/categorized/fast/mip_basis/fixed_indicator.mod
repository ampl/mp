# Indicator constraint inactive at the optimum (b=0).
# Gurobi's fixed model drops it by default (FixVarsInIndicators=0),
# so y stays basic in the fixed LP; with FixVarsInIndicators=1
# y is fixed to its solution value and becomes nonbasic.
var b binary;
var y >= 0, <= 10;
var z >= 0;
minimize cost: 3*b + y + 2*z;
subject to demand: y + z >= 2;
subject to logic: b = 1 ==> y <= 1;
option presolve 0;
