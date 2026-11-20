var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (-1)*x[2] + (0)*x[3] + (2)*(if abs(x[1] - x[3]) < 0 then 1 else 0);
s.t. c1: not (max(x[1], x[3]) > 0);
s.t. c2: (if x[2] > 0 then x[2] else 0) + (if max(x[1], x[3]) > 0 then x[2] else 0) <= 0;
s.t. c3: (if x[2] > 0 then x[3] else x[3]) > 2;
