var x{1..6} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (-2)*x[2] + (4)*x[3] + (6)*x[4] + (-1)*x[5] + (6)*x[6];
s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 3;
s.t. c2: (if x[3] == 4 then x[4] else x[2]) != 4;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 0;
