var x{1..4} integer >= -2 <= 1;
maximize o: (-2)*x[1] + (-4)*x[2] + (6)*x[3] + (3)*x[4];
s.t. c1: alldiff(x[3], x[4], x[1]);
s.t. c2: max(x[3], x[1]) - min(x[4], x[1]) != -1;
s.t. c3: numberof -1 in (x[1], x[2], x[3], x[4]) != 1;
