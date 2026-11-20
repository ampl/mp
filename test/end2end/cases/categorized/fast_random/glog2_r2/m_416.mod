var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 2;
s.t. c2: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
