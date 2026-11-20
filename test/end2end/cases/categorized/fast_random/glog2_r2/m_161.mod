var x{1..4} integer >= -2 <= 0;
maximize o: (0)*x[1] + (1)*x[2] + (-2)*x[3] + (-4)*x[4];
s.t. c1: ((if x[1] < -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
