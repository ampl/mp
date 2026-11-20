var x{1..5} integer >= 0 <= 2;
maximize o: (2)*x[1] + (-1)*x[2] + (-2)*x[3] + (2)*x[4] + (6)*x[5];
s.t. c1: abs(x[2] - x[4]) < 1;
