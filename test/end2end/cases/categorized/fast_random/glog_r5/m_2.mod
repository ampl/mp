var x{1..5} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-4)*x[2] + (-2)*x[3] + (-4)*x[4] + (-3)*x[5];
s.t. c1: max(x[2], x[3]) - min(x[2], x[5]) <= 1;
s.t. c2: (x[4] > 2) or (x[4] != 0) or (x[5] < 1);
s.t. c3: abs(x[4] - x[3]) == 0;
