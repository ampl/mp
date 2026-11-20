var x{1..5} integer >= -3 <= 1;
maximize o: (0)*x[1] + (6)*x[2] + (5)*x[3] + (-3)*x[4] + (-1)*x[5];
s.t. c1: ((if x[3] == 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
