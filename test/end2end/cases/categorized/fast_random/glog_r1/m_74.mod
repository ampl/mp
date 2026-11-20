var x{1..4} integer >= -2 <= 0;
minimize o: (-1)*x[1] + (-3)*x[2] + (0)*x[3] + (5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -2;
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4]) <= 3;
