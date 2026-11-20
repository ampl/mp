var x{1..4} integer >= 0 <= 4;
maximize o: (4)*x[1] + (3)*x[2] + (-2)*x[3] + (3)*x[4];
s.t. c1: atleast 3 {i in 1..4} (x[i] >= 4);
s.t. c2: ((if x[3] == 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
