var x{1..4} integer >= -2 <= 0;
maximize o: (1)*x[1] + (-4)*x[2] + (-3)*x[3] + (4)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: (x[4] >= -2) ==> (x[3] > -2);
s.t. c3: ((if x[3] != -1 then 1 else 0) + (if x[1] >= -1 then 1 else 0)) >= 1;
