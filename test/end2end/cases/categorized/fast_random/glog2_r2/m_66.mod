var x{1..3} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: max(x[3], x[1]) - min(x[2], x[1]) > -2;
s.t. c2: ((if x[2] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
