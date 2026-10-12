var b binary;
var y >= 0;
minimize cost: b + y;
subject to logic: b = 0 ==> y >= 2;
option presolve 0;
