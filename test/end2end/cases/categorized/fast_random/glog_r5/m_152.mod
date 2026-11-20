var x{1..5} integer >= 0 <= 4;
maximize o: (3)*x[1] + (2)*x[2] + (-1)*x[3] + (0)*x[4] + (6)*x[5];
s.t. c1: abs(x[3] - x[5]) == 1;
