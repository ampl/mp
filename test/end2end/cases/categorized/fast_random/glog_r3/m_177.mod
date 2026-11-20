var x{1..5} integer >= 0 <= 3;
maximize o: (3)*x[1] + (2)*x[2] + (-4)*x[3] + (5)*x[4] + (0)*x[5];
s.t. c1: not ((x[4] <= 3) and (x[2] >= 2));
s.t. c2: abs(x[5] - x[2]) >= 1;
