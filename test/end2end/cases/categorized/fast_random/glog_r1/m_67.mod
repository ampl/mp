var x{1..4} integer >= 0 <= 3;
maximize o: (1)*x[1] + (1)*x[2] + (0)*x[3] + (6)*x[4];
s.t. c1: (x[4] <= 0) ==> (x[4] == 0);
s.t. c2: ((if x[3] != 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
