var x{1..4} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (5)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: ((if x[2] > -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
