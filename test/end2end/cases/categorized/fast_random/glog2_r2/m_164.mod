var x{1..5} integer >= -2 <= 3;
minimize o: (-1)*x[1] + (-3)*x[2] + (4)*x[3] + (1)*x[4] + (4)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 3;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) > 1;
