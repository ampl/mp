var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (-3)*x[2] + (6)*x[3];
s.t. c1: max(x[2], x[2]) - min(x[1], x[2]) <= 2;
s.t. c2: abs(x[3] - x[1]) < 2;
s.t. c3: x[2] == 0;
