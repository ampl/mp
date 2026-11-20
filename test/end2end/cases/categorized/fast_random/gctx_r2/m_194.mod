var x{1..5} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (-4)*x[2] + (-3)*x[3] + (0)*x[4] + (5)*x[5] + (-2)*(if x[2] + x[3] <= 2 then 1 else 0);
s.t. c1: (x[2] + x[3] <= 2) or (max(x[1], x[3]) > 2);
s.t. c2: not ((max(x[1], x[3]) > 2) and (max(x[1], x[3]) > 2));
