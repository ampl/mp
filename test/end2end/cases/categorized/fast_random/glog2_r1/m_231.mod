var x{1..4} integer >= -2 <= 0;
maximize o: (1)*x[1] + (-3)*x[2] + (2)*x[3] + (0)*x[4];
s.t. c1: abs(x[3] - x[4]) > 1;
