var x{1..5} integer >= -3 <= 1;
maximize o: (2)*x[1] + (6)*x[2] + (6)*x[3] + (-2)*x[4] + (1)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c3: max(x[1], x[4]) - min(x[3], x[1]) != -3;
s.t. c4: abs(x[4] - x[4]) < 0;
