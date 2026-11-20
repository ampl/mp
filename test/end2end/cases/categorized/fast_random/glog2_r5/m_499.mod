var x{1..5} integer >= 0 <= 3;
maximize o: (1)*x[1] + (1)*x[2] + (6)*x[3] + (1)*x[4] + (-5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
