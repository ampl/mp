var x{1..5} integer >= 0 <= 4;
maximize o: (4)*x[1] + (6)*x[2] + (4)*x[3] + (-2)*x[4] + (6)*x[5];
s.t. c1: abs(x[3] - x[1]) == 4;
