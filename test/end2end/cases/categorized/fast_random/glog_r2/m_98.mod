var x{1..4} integer >= 0 <= 2;
maximize o: (4)*x[1] + (1)*x[2] + (2)*x[3] + (-4)*x[4];
s.t. c1: ((if x[3] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..4} x[i] == 1;
