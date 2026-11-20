var x{1..4} integer >= -2 <= 3;
maximize o: (3)*x[1] + (-1)*x[2] + (0)*x[3] + (2)*x[4];
s.t. c1: (x[3] == 2) <==> (x[1] + x[2] <= 0);
s.t. c2: ((if x[1] <= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
