var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-2)*x[2] + (-2)*x[3] + (4)*x[4] + (5)*x[5] + (0)*(if x[5] == 2 then 1 else 0);
s.t. c1: (min(x[1], x[1]) < 1) ==> ((x[5] == 2) or (x[3] > 1));
s.t. c2: (x[3] > 1) and ((min(x[1], x[1]) < 1) or not (min(x[1], x[1]) < 1));
s.t. c3: (if x[3] > 1 then x[4] else 0) + (if min(x[1], x[1]) < 1 then x[1] else 0) <= 0;
