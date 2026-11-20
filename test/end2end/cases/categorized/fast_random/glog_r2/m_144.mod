var x{1..6} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (-2)*x[2] + (4)*x[3] + (4)*x[4] + (5)*x[5] + (6)*x[6];
s.t. c1: max(x[4], x[4]) - min(x[6], x[2]) > 2;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
