var x{1..3} integer >= -2 <= 1;
maximize o: (0)*x[1] + (-4)*x[2] + (6)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: not ((x[3] >= 1) and (x[1] <= -1));
s.t. c3: ((if x[2] > 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
