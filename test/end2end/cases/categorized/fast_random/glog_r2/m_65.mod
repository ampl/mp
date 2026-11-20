var x{1..4} integer >= -2 <= 1;
maximize o: (1)*x[1] + (-4)*x[2] + (6)*x[3] + (6)*x[4];
s.t. c1: abs(x[3] - x[1]) < 2;
