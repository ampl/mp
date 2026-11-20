var x{1..3} integer >= 0 <= 2;
maximize o: (4)*x[1] + (5)*x[2] + (6)*x[3];
s.t. c1: max(x[1], x[3]) - min(x[3], x[3]) != 2;
s.t. c2: (x[3] != 2) ==> (x[3] > 0);
s.t. c3: (x[3] != 0) or (x[3] <= 1) or (x[1] >= 0);
