var x{1..5} integer >= -2 <= 0;
maximize o: (4)*x[1] + (-3)*x[2] + (4)*x[3] + (0)*x[4] + (1)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: not ((x[4] <= -3) and (x[5] > 0));
