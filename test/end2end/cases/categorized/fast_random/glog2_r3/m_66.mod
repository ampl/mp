var x{1..3} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: not ((x[2] <= 2) and (x[1] > 1));
s.t. c2: ((if x[1] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: max(x[1], x[1]) - min(x[2], x[2]) == 0;
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
