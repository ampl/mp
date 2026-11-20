var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (1)*x[2] + (-1)*x[3] + (-2)*x[4];
s.t. c1: (min(x[3], x[4]) >= 2) or (max(x[4], x[2]) > 3);
s.t. c2: (min(x[3], x[2]) <= 0) and ((min(x[3], x[4]) >= 2) or not (max(x[4], x[2]) > 3));
s.t. c3: (min(x[3], x[2]) <= 0) <==> (max(x[4], x[2]) > 3);
