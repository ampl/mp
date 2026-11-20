var x{1..6} integer >= 0 <= 3;
maximize o: (5)*x[1] + (1)*x[2] + (-2)*x[3] + (0)*x[4] + (4)*x[5] + (3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[5] > 2) and (x[2] < 0) ==> (x[6] >= 2);
s.t. c3: (x[1] >= 1) or (x[5] <= 3) or (x[1] >= 2);
