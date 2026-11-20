var x{1..5} integer >= -2 <= 0;
maximize o: (3)*x[1] + (-3)*x[2] + (-4)*x[3] + (-3)*x[4] + (4)*x[5];
s.t. c1: ((if x[5] != -2 then 1 else 0) + (if x[5] >= -2 then 1 else 0)) >= 1;
