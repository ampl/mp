var x{1..6} integer >= -2 <= 1;
maximize o: (4)*x[1] + (3)*x[2] + (5)*x[3] + (-3)*x[4] + (-4)*x[5] + (3)*x[6];
s.t. c1: ((if x[1] <= -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
