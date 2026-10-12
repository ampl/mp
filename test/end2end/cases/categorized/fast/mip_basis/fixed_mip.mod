# The fixed LP has demand dual 3; the LP relaxation has dual 3.5.
var b binary := 1;
var k integer >= 0, <= 4 := 2;
var y >= 0, <= 2.5 := 0.25;
minimize cost: b + 7*k + 3*y;
subject to demand: b + 2*k + y >= 5.25;
option presolve 0;
