var x{1..5} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: ((if x[2] != 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
