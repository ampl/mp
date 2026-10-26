# -------------------------------------------------------------
# prod_int_01.mod
# Products of integer variables: bilinear terms, a product
# of sums, integer powers, a quadratic objective, and
# the same product x*y in the objective and a constraint
# (CP solvers share one product constraint).
# -------------------------------------------------------------

var x integer in -5..9;
var y integer in -3..7;
var z integer in 0..6;
var w integer in 0..6;
maximize Obj: x*y - z^3 + (z+w)*(x-y) - 2*x*y;          # quadratic objective, repeated x*y
s.t. P1: x*y <= 20;
s.t. P2: (x+z)*(y-w) >= -15;                             # product of sums
s.t. C3: z^3 + w^2 <= 60;                                # integer powers
