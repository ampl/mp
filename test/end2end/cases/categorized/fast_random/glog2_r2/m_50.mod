var x{1..4} integer >= 0 <= 2;
maximize o: (6)*x[1] + (5)*x[2] + (4)*x[3] + (2)*x[4];
s.t. c1: abs(x[1] - x[3]) != 1;
