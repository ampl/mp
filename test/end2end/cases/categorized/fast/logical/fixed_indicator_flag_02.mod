# -------------------------------------------------------------
# fixed_indicator_flag_02.mod
# Indicator flag (x[4] < 0) fixed after the indicators of c1
# were created (cvt:pre:indfixedflag).
# -------------------------------------------------------------

var x{1..6} integer >= -2 <= 2;
maximize o: (3)*x[1] + (5)*x[2] + (-5)*x[3] + (-1)*x[4] + (6)*x[5] + (2)*x[6];
s.t. c1: (if x[4] < 0 then x[5] else x[6]) <= -1;
s.t. c2: (x[4] <= 2) ==> (x[4] >= 0) else (x[5] < -2);
