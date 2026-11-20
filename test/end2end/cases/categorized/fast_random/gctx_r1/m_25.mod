var x{1..4} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-4)*x[2] + (-3)*x[3] + (-2)*x[4];
s.t. c1: not (min(x[4], x[4]) == 1);
s.t. c2: (if min(x[4], x[4]) == 1 then x[3] else x[4]) > 0;
s.t. c3: (min(x[4], x[4]) == 1) and ((x[3] <= 2) or not (min(x[4], x[2]) <= 1));
s.t. c4: (min(x[4], x[4]) == 1) ==> (x[3] <= 2);
