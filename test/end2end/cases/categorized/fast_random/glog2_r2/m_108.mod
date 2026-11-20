var x{1..6} integer >= 0 <= 1;
maximize o: (3)*x[1] + (1)*x[2] + (-4)*x[3] + (5)*x[4] + (2)*x[5] + (0)*x[6];
s.t. c1: ((if x[2] < 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
