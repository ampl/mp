# Unique conflict; removing any row makes the model feasible.
# Two-variable rows prevent conversion from resolving it by bound propagation.
var x;
var y;
var b binary;
var z;
var c binary;
subject to Cap: x + y <= 0;
# The first indicator is excluded from the conflict. Mixing its index with
# linear-row indices would incorrectly overwrite Cap's membership.
subject to Irrelevant: c = 1 ==> z >= 0;
subject to WhenZero: b = 0 ==> x + y >= 1;
subject to WhenOne: b = 1 ==> x + y >= 2;
option presolve 0;
