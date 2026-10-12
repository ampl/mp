var b binary;
var x;
minimize cost: b;
subject to lower: x >= 2*b + 1;
subject to upper: x <= 2*b;
option presolve 0;
