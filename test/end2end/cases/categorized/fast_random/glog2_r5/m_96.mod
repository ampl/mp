var x{1..4} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: max(x[1], x[2]) - min(x[1], x[1]) < -1;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: ((x[3] != 1) ==> (x[2] == 0)) and ((x[2] == 0) ==> (x[3] >= 1));
