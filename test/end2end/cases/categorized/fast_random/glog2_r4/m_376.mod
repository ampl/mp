var x{1..4} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: abs(x[4] - x[2]) != 2;
s.t. c2: numberof 3 in (x[1], x[2], x[3], x[4]) < 1;
s.t. c3: (x[4] > 0) or (x[2] >= 1) or (x[1] == 2);
