var x{1..6} integer >= 0 <= 4;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (if x[2] <= 3 then x[1] else x[3]) > 2;
s.t. c3: ((if x[1] != 0 then 1 else 0) + (if x[6] >= 0 then 1 else 0)) >= 1;
