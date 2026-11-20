var x{1..5} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (4)*x[2] + (0)*x[3] + (3)*x[4] + (2)*x[5];
s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c2: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) < 3;
