# Continuous variables in an SOS1 constraint, primal start values only.
# No integer variables, but solvers with native SOS solve a MIP and take
# the values as a MIP start: no warning about unused start values.
# Optimum 30: x[3] = 10.
suffix sosno IN;
suffix ref IN;

var x {1..3} >= 0 <= 10 := 1;

maximize o: x[1] + 2*x[2] + 3*x[3];

s.t. c: sum {i in 1..3} x[i] <= 15;

let {i in 1..3} x[i].sosno := 1;    # one SOS1 set
let {i in 1..3} x[i].ref := i;
