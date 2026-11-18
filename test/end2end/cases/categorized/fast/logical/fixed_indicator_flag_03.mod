# -------------------------------------------------------------
# fixed_indicator_flag_03.mod
# Indicator flags fixed after the indicators of c3
# were created (cvt:pre:indfixedflag).
# -------------------------------------------------------------

var x{1..4} integer >= -1 <= 2;
minimize o: (3)*x[1] + (-4)*x[2] + (4)*x[3] + (-2)*x[4];
s.t. c1: (abs(x[2] - x[3]) == 0) <==> (x[3] != -1);
s.t. c2: (abs(x[2] - x[3]) == 0) ==> (x[3] != -1) else (abs(x[2] - x[3]) == 0);
s.t. c3: (if x[3] + x[1] == 1 then x[3] else 0) + (if abs(x[2] - x[3]) == 0 then x[3] else 0) <= 1;
