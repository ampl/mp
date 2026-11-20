var x{1..6} integer >= -2 <= 2;
minimize o: (0)*x[1] + (-2)*x[2] + (1)*x[3] + (-5)*x[4] + (2)*x[5] + (-3)*x[6];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 0;
s.t. c2: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 3;
