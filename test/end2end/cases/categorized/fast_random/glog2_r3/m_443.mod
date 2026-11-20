var x{1..3} integer >= 0 <= 4;
maximize o: (5)*x[1] + (1)*x[2] + (6)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: abs(x[2] - x[2]) > 0;
s.t. c3: max(x[2], x[3]) - min(x[1], x[1]) > 4;
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
s.t. c5: x[2] == 1;
s.t. c6: x[1] == 4;
