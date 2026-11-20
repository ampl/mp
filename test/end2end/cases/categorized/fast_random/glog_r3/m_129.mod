var x{1..6} integer >= 0 <= 4;
maximize o: (-5)*x[1] + (6)*x[2] + (2)*x[3] + (-3)*x[4] + (2)*x[5] + (-4)*x[6];
s.t. c1: (if x[6] >= 3 then x[4] else x[6]) <= 0;
