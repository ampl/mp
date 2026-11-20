var x{1..6} integer >= -2 <= 1;
minimize o: (6)*x[1] + (6)*x[2] + (-2)*x[3] + (-1)*x[4] + (4)*x[5] + (5)*x[6];
s.t. c1: (if x[6] > 0 then x[6] else x[5]) == -1;
