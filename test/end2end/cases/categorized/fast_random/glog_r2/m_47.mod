var x{1..5} integer >= 0 <= 2;
minimize o: (2)*x[1] + (-2)*x[2] + (-2)*x[3] + (-3)*x[4] + (-4)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 1;
