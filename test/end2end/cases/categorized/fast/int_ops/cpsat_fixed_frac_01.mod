# -------------------------------------------------------------
# cpsat_fixed_frac_01.mod
# A continuous variable fixed at a fractional value: CP-SAT has
# only integer variables; the driver passes it as a constant.
# -------------------------------------------------------------
var x >= 2.5 <= 2.5;
var y integer >= 0 <= 5;
minimize o: x + 3*y;
s.t. c: y + x >= 4;
