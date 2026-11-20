var x{1..4} integer >= -2 <= 2;
maximize o: (1)*x[1] + (3)*x[2] + (-4)*x[3] + (0)*x[4];
s.t. c1: ((if x[1] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: not ((x[4] != 2) and (x[4] > 1));
s.t. c3: (x[1] > 0) ==> (x[4] <= 0);
