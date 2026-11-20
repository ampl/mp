var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (-4)*x[2] + (5)*x[3] + (-2)*x[4] + (-1)*(if x[4] + x[1] == 0 then 1 else 0);
s.t. c1: not ((max(x[1], x[4]) == 1) and (x[4] + x[1] == 0));
s.t. c2: (max(x[1], x[4]) == 1) ==> ((max(x[2], x[3]) < 3) or (max(x[1], x[4]) == 1));
