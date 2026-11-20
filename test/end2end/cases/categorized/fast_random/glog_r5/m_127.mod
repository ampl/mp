var x{1..3} integer >= -2 <= 2;
maximize o: (0)*x[1] + (6)*x[2] + (-5)*x[3];
s.t. c1: max(x[3], x[1]) - min(x[1], x[2]) != -2;
s.t. c2: ((if x[2] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[2] >= -2) ==> (x[1] <= 0);
