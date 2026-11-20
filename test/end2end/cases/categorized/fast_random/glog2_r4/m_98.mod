var x{1..4} integer >= 0 <= 1;
maximize o: (3)*x[1] + (3)*x[2] + (3)*x[3] + (4)*x[4];
s.t. c1: ((if x[2] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: max(x[4], x[3]) - min(x[3], x[4]) > 1;
s.t. c3: atleast 2 {i in 1..4} (x[i] >= 0);
