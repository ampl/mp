var x{1..6} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-4)*x[2] + (4)*x[3] + (-1)*x[4] + (1)*x[5] + (6)*x[6];
s.t. c1: (x[4] == 0) and (x[4] != 1) ==> (x[1] >= 0);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: max(x[6], x[2]) - min(x[3], x[4]) <= 2;
