var x{1..5} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
