# -------------------------------------------------------------
# pow_int_01.mod
# Integer powers and products of integer variables
# -------------------------------------------------------------

var a integer in -4..4;
var b integer in 0..5;
var c integer in -3..3;

maximize Obj: a^3 - b^2 + a*b - c^4 + 3*c;
s.t. Cubes: a^2 + b^3 <= 40;
s.t. Prod3: a*b*c >= -6;             # product of 3 variables
