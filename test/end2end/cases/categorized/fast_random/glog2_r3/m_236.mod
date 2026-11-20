var x{1..4} integer >= -2 <= 1;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (if x[2] <= -1 then x[4] else x[2]) > 0;
s.t. c2: (x[1] <= -2) ==> (x[1] <= -1) else (x[2] <= -2);
s.t. c3: ((if x[4] > -1 then 1 else 0) + (if x[1] >= -1 then 1 else 0)) >= 1;
