var x{1..4} integer >= 0 <= 2;
minimize o: (2)*x[1] + (-5)*x[2] + (4)*x[3] + (2)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) == 3;
s.t. c2: exists {i in 1..4} x[i] == 0;
