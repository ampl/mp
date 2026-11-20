var x{1..4} integer >= 0 <= 2;
maximize o: (6)*x[1] + (4)*x[2] + (2)*x[3] + (-3)*x[4];
s.t. c1: (x[4] + x[3] <= 0) ==> ((x[4] + x[3] <= 0) or (x[4] + x[3] <= 0));
s.t. c2: (if min(x[1], x[3]) == 1 then x[2] else x[3]) < 2;
