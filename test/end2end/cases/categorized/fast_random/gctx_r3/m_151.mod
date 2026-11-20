var x{1..4} integer >= 0 <= 3;
maximize o: (0)*x[1] + (-2)*x[2] + (-3)*x[3] + (-2)*x[4] + (-2)*(if max(x[4], x[4]) <= 2 then 1 else 0);
s.t. c1: (x[2] >= 1) ==> ((x[2] >= 1) or (x[4] > 0));
s.t. c2: (x[4] > 0) and ((x[2] >= 1) or not (x[2] >= 1));
s.t. c3: (x[4] > 0) ==> (max(x[4], x[4]) <= 2);
s.t. c4: not (max(x[4], x[4]) <= 2);
