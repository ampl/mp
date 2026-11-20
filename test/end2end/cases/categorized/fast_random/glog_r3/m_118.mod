var x{1..4} integer >= -2 <= 2;
maximize o: (5)*x[1] + (6)*x[2] + (2)*x[3] + (6)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
