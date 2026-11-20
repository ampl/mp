var x{1..4} integer >= 0 <= 2;
maximize o: (3)*x[1] + (2)*x[2] + (1)*x[3] + (5)*x[4];
s.t. c1: not ((min(x[3], x[3]) == 2) and (x[2] + x[1] >= 0));
s.t. c2: (min(x[3], x[3]) == 2) or (x[2] + x[1] >= 0);
s.t. c3: (min(x[3], x[3]) == 2) or (x[2] + x[1] >= 0);
