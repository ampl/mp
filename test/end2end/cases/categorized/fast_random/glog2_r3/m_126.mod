var x{1..4} integer >= 0 <= 4;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: abs(x[1] - x[3]) >= 1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: abs(x[3] - x[2]) >= 2;
s.t. c4: (x[3] > 4) or (x[2] >= 1) or (x[2] >= 3);
