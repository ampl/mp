var x{1..4} integer >= 0 <= 3;
maximize o: (6)*x[1] + (5)*x[2] + (3)*x[3] + (2)*x[4] + (0)*(if max(x[3], x[4]) != 1 then 1 else 0);
s.t. c1: (max(x[3], x[4]) != 1) ==> ((x[2] > 0) or (max(x[3], x[4]) != 1));
s.t. c2: (x[2] > 0) and ((x[1] + x[4] < 0) or not (x[2] > 0));
