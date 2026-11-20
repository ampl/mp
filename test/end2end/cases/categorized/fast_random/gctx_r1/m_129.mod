var x{1..3} integer >= 0 <= 2;
maximize o: (2)*x[1] + (-3)*x[2] + (3)*x[3] + (-6)*(if max(x[1], x[1]) <= 2 then 1 else 0);
s.t. c1: (if max(x[1], x[1]) <= 2 then x[3] else 0) + (if max(x[1], x[1]) <= 2 then x[2] else 0) >= 2;
s.t. c2: not (max(x[1], x[1]) <= 2);
