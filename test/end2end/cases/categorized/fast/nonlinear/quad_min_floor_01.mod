# A constraint with a product and a term without derivatives
# (min of a constant and floor): not quadratic.  degree_ASL used to call
# it quadratic, the quadratic decomposition then failed, and Knitro
# stopped.
var x := 1.3; var y := 2.1; var z := 0.7;

minimize o: x^2 + y^2;

s.t. c: x*y + min(1, floor(z)) <= 30;
