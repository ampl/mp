var x{1..6} integer >= 0 <= 2;
maximize o: (3)*x[1] + (4)*x[2] + (2)*x[3] + (-5)*x[4] + (6)*x[5] + (-5)*x[6];
s.t. c1: (if x[1] > 0 then x[5] else x[2]) != 1;
