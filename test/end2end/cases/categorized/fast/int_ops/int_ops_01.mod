# -------------------------------------------------------------
# int_ops_01.mod
# Integer division and remainder (AMPL div, mod), also with
# variable and negative divisors, and rounding of
# integer expressions, together in one integer model.
# -------------------------------------------------------------

var x integer in -20..20;
var y integer in -5..5;          # divisor: 0 is excluded by div
var z integer in -20..20;
var w integer in 0..9;

s.t. Quot:  x div y >= 3;                 # variable divisor
s.t. Rem:   x mod y = -2;                 # sign of x
s.t. Rem4:  z mod 4 = 1;
s.t. Mixed: z div 3 + x div 4 - w <= 1;   # constant divisors
s.t. RoundInt: round(w/1) + trunc(z) <= 12;

maximize Obj: 10*x - y + 2*z + (x + z) mod 7 + w;
