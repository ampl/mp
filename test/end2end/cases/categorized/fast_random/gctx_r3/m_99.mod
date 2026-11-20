var x{1..4} integer >= 0 <= 3;
maximize o: (2)*x[1] + (-5)*x[2] + (6)*x[3] + (-1)*x[4] + (2)*(if x[3] >= 2 then 1 else 0);
s.t. c1: (x[2] + x[4] == 1) or (max(x[2], x[2]) > 1);
s.t. c2: not ((max(x[2], x[2]) > 1) and (x[2] + x[4] == 1));
s.t. c3: not (max(x[2], x[2]) > 1);
