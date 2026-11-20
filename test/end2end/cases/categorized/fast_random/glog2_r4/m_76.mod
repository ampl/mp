var x{1..6} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: (x[5] != 1) and (x[1] != 2) ==> (x[4] >= 1);
s.t. c2: ((if x[2] <= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
