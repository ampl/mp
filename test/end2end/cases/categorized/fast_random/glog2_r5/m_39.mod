var x{1..6} integer >= -2 <= 0;
minimize o: (1)*x[1] + (-4)*x[2] + (-3)*x[3] + (2)*x[4] + (-5)*x[5] + (-5)*x[6];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
s.t. c2: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 3;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
s.t. c4: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 3;
