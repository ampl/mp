var x{1..3} integer >= -1 <= 2;
maximize o: (-2)*x[1] + (-2)*x[2] + (2)*x[3] + (3)*(if max(x[2], x[3]) > -1 then 1 else 0);
s.t. c1: (x[1] > 1) and ((x[1] > 1) or not (max(x[2], x[3]) > -1));
s.t. c2: (if min(x[1], x[3]) > 2 then x[2] else x[2]) <= 2;
