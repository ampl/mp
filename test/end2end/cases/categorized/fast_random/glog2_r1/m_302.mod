var x{1..6} integer >= 0 <= 3;
maximize o: (4)*x[1] + (0)*x[2] + (-4)*x[3] + (-1)*x[4] + (1)*x[5] + (4)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
