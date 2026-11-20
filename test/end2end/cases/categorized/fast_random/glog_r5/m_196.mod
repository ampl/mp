var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (6)*x[2] + (1)*x[3];
s.t. c1: (x[3] < 2) or (x[2] == 2) or (x[2] > 2);
s.t. c2: abs(x[3] - x[3]) != 1;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: exists {i in 1..3} x[i] == 1;
