var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (-3)*x[2] + (5)*x[3] + (1)*x[4];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) != 2;
s.t. c2: ((if x[4] == 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
