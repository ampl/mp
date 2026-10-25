# Chained PL terms: the argument of a PL term is a linear
# expression of another PL term's result. All variables are
# unbounded; the arguments are fixed by constraints.
# Needs AMPL presolve off and pl_linearize=0 to reach the solver.
# Without result bounds of PL terms, the PL -> SOS2 -> MIP
# reformulation extended the end segments to +-1e6, and the
# solver's feasibility tolerance gave objective 0.

var x1;  s.t. pin_x1: x1 = -1.5;
var x2;  s.t. pin_x2: x2 = 0;
var x3;  s.t. pin_x3: x3 = 0.5;

var y1 = -2*<<-1; -1,-0.5>> x1 - x2 + 1;               # -1
var y2 = <<-4,2,4; 3,0.5,3,-1>> y1 + 0*x3;             # -0.5

var t;
maximize o: t;                                         # -0.75
s.t. ct: t <= -<<-2,1,4; -0.5,0.5,0.5,0.5>> y2 + y1;
