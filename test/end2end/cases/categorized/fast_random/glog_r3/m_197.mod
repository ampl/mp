var x{1..3} integer >= 0 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: numberof 1 in (x[1], x[2], x[3]) == 3;
s.t. c2: numberof 2 in (x[1], x[2], x[3]) <= 2;
