var x{1..4} integer >= -1 <= 2;
maximize o: (5)*x[1] + (2)*x[2] + (-2)*x[3] + (-1)*x[4] + (-2)*(if min(x[2], x[1]) != 0 then 1 else 0);
s.t. c1: (min(x[2], x[1]) != 0) ==> (abs(x[1] - x[1]) < 3);
s.t. c2: (max(x[1], x[3]) < 2) or (x[2] >= 2);
