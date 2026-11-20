var x{1..5} integer >= -2 <= 1;
minimize o: (4)*x[1] + (2)*x[2] + (2)*x[3] + (4)*x[4] + (-4)*x[5];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
s.t. c2: abs(x[1] - x[2]) >= 1;
