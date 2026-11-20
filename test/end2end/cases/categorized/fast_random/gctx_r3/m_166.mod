var x{1..4} integer >= -1 <= 1;
maximize o: (-2)*x[1] + (-1)*x[2] + (-2)*x[3] + (4)*x[4];
s.t. c1: not ((min(x[3], x[3]) <= -1) and (min(x[2], x[1]) == 0));
s.t. c2: (min(x[3], x[3]) <= -1) or (min(x[3], x[3]) <= -1);
s.t. c3: (min(x[2], x[1]) == 0) and ((min(x[2], x[1]) == 0) or not (min(x[3], x[3]) <= -1));
