var x{1..4} integer >= 0 <= 4;
maximize o: (2)*x[1] + (1)*x[2] + (2)*x[3] + (5)*x[4];
s.t. c1: abs(x[3] - x[3]) != 3;
