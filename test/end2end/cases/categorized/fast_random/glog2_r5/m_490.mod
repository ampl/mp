var x{1..4} integer >= -2 <= 0;
maximize o: (2)*x[1] + (-2)*x[2] + (0)*x[3] + (3)*x[4];
s.t. c1: ((if x[2] > -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
