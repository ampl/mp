#####################################################
# A quadratic and a linear objective with equal
# priorities: native multi-objective handling must not
# drop the quadratic one silently (COPT did).
#####################################################

var x >= 0 <= 10;
var y >= 0 <= 10;

minimize Obj1QP: 2*x^2 + 0.5*y^2 - 4*x - 2*y;
minimize Obj2Lin: x + y;
