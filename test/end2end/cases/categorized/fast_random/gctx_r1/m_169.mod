var x{1..3} integer >= 0 <= 3;
maximize o: (3)*x[1] + (4)*x[2] + (6)*x[3];
s.t. c1: (min(x[2], x[2]) < 3) ==> (x[2] >= 1);
s.t. c2: (if x[3] <= 3 then x[1] else 0) + (if min(x[2], x[2]) < 3 then x[3] else 0) > 3;
