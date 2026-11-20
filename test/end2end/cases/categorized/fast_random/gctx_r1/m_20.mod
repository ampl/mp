var x{1..5} integer >= -1 <= 2;
maximize o: (3)*x[1] + (-4)*x[2] + (0)*x[3] + (-3)*x[4] + (6)*x[5];
s.t. c1: not ((x[4] < -1) and (x[5] != -1));
s.t. c2: not (max(x[1], x[1]) > 0);
