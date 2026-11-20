var x{1..4} integer >= 0 <= 4;
maximize o: (5)*x[1] + (-4)*x[2] + (2)*x[3] + (4)*x[4];
s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4]) < 0;
s.t. c2: abs(x[2] - x[1]) != 3;
