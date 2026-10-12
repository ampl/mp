suffix sosno integer IN;
var v {1..2} >= 0, <= 1;
var b binary;
minimize cost: b + v[1] + v[2];
subject to demand: v[1] + v[2] >= 1;
subject to link: b >= v[1];
let v[1].sosno := 1;
let v[2].sosno := 1;
option presolve 0;
