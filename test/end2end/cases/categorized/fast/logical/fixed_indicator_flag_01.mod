# -------------------------------------------------------------
# fixed_indicator_flag_01.mod
# The flag of the indicators of c1 (x[4] <= 1) is fixed by c2
# only after the indicators were created (cvt:pre:indfixedflag).
# Some solvers reject a constant (continuous) indicator flag.
# -------------------------------------------------------------

var x{1..4} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-3)*x[2] + (6)*x[3] + (4)*x[4];
s.t. c1: (if x[4] <= 1 then x[2] else x[4]) <= 0;
s.t. c2: (x[4] < 2) ==> (max(x[2], x[4]) != 0) else (x[4] <= 1);
s.t. c3: (min(x[1], x[2]) >= 0) <==> (max(x[2], x[4]) != 0);
