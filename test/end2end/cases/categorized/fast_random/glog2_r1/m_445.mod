var x{1..3} integer >= 0 <= 3;
maximize o: (0)*x[1] + (2)*x[2] + (3)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] > 1) or (x[1] == 1) or (x[3] == 1);
s.t. c3: not ((x[1] >= 2) and (x[3] >= 3));
s.t. c4: exists {i in 1..3} x[i] == 0;
s.t. c5: max(x[2], x[3]) - min(x[1], x[1]) != 3;
s.t. c6: (x[1] > 0) or (x[1] > 2) or (x[1] == 0);
